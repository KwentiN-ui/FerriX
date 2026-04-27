//! Static stress analysis step.
//!
//! This module implements the solver for static FEA problems, supporting both
//! linear and non-linear (NLGEOM) analyses using the Newton-Raphson method.

use crate::solver::assembler::Assembler;
use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::NodeId;
use crate::solver::increment::IncrementData;
use crate::solver::io::writer::ResultWriter;
use crate::solver::preconditioner::DiagonalPreconditioner;
use crate::solver::project::Project;
use crate::solver::results::{FieldType, IncResult, NodalResult};
use crate::solver::solvers::SolverType;
use crate::solver::solvers::{Solver, direct::DirectSolver, iterative::IterativeSolver};
use crate::solver::state::SolutionState;
use crate::solver::step::boundary_conds::{BoundaryCondition, Load};
use crate::solver::time::SolverTime;
use sprs::CsMat;
use std::collections::HashMap;

/// Configuration and state for a static analysis step.
#[derive(Debug, Clone)]
pub struct StaticStep {
    /// The linear solver strategy to use.
    pub solver: SolverType,
    /// Time incrementation settings.
    pub increment_data: IncrementData,
    /// Concentrated loads applied during this step.
    pub loads: Vec<Load>,
    /// Boundary conditions applied during this step.
    pub bcs: Vec<BoundaryCondition>,
    /// If true, performs a non-linear analysis considering geometric non-linearity.
    pub nlgeom: bool,
}

impl StaticStep {
    /// Executes the static analysis by iterating through increments.
    ///
    /// For non-linear analyses, it employs the Newton-Raphson iteration loop
    /// within each increment to find the equilibrium state.
    ///
    /// # Errors
    /// Returns `FerrixError::ConvergenceError` if the solver fails to converge
    /// within the allowed iterations or if the minimum time increment is reached.
    pub fn solve(
        &self,
        step_id: usize,
        project: &Project,
        solution_state: &mut SolutionState,
        writer: &dyn ResultWriter,
        timer: &mut SolverTime,
    ) -> Result<()> {
        if self.nlgeom {
            println!("--- Step {step_id}: StaticStep (Non-Linear) ---");
        } else {
            println!("--- Step {step_id}: StaticStep (Linear) ---");
        }

        let mut current_time = 0.0;
        let mut n_inc = 0;
        let mut increment_time = self.increment_data.initial_time_increment;

        // The displacement at the start of the step
        let mut u_step_start = solution_state.displacements.clone();

        while current_time < self.increment_data.time_period {
            n_inc += 1;
            if n_inc > self.increment_data.max_iterations {
                return Err(FerrixError::ConvergenceError(format!(
                    "Maximum increment count of {} exceeded.",
                    self.increment_data.max_iterations,
                )));
            }

            if current_time + increment_time > self.increment_data.time_period {
                increment_time = self.increment_data.time_period - current_time;
            }

            timer.new_increment(increment_time);

            println!(
                "Increment {n_inc} | Step Time: {:.4e} | Increment Size: {increment_time:.4e}",
                timer.local_time()
            );

            // --- Newton-Raphson Loop ---
            let mut converged = false;
            let mut u_curr = u_step_start.clone(); // Start iteration with displacement from end of previous increment
            let num_dofs = u_curr.len();

            for iter in 0..15 {
                // Max NR iterations
                // 1. Calculate External Forces F_ext
                let mut f_ext = vec![0.0; num_dofs];
                for load in &self.loads {
                    if let Some(idx) = project.mesh.get_index_for_node_id(load.node_id()) {
                        let global_dof = idx * 3 + load.dof();
                        f_ext[global_dof] += load.value(timer, step_id);
                    }
                }

                // 2. Calculate Internal Forces F_int
                // If NLGEOM is active, internal forces are integrated over the current configuration.
                // Otherwise (Linear), we integrate over the initial configuration.
                let u_conf = if self.nlgeom {
                    Some(u_curr.as_slice())
                } else {
                    None
                };
                let f_int = Assembler::assemble_internal_force(project, &u_curr, u_conf)?;

                // 3. Calculate Residual R = F_ext - F_int
                let mut residual = vec![0.0; num_dofs];
                for i in 0..num_dofs {
                    residual[i] = f_ext[i] - f_int[i];
                }

                // 4. Check Convergence
                // For the convergence check, we only consider DOFs that do not have a boundary condition applied.
                // Boundary conditions generate reaction forces in F_int, which do not exist in F_ext,
                // so the residual R = F_ext - F_int will never be zero at those nodes.
                let mut res_for_check = residual.clone();
                let mut f_ext_for_check = f_ext.clone();
                for bc in &self.bcs {
                    if let Some(idx) = project.mesh.get_index_for_node_id(bc.node_id()) {
                        let global_dof = idx * 3 + bc.dof();
                        if global_dof < num_dofs {
                            res_for_check[global_dof] = 0.0;
                            f_ext_for_check[global_dof] = 0.0;
                        }
                    }
                }

                let res_norm: f64 = res_for_check.iter().map(|&r| r * r).sum::<f64>().sqrt();
                let f_ext_norm: f64 = f_ext_for_check.iter().map(|&f| f * f).sum::<f64>().sqrt();

                let rel_res = if f_ext_norm > 1e-6 {
                    res_norm / f_ext_norm
                } else {
                    res_norm
                };

                println!("  Iteration {iter}: |R| = {res_norm:.3e}, |R|/|F_ext| = {rel_res:.3e}");

                if rel_res < 1e-3 {
                    // Convergence tolerance
                    converged = true;
                    break;
                }

                // 5. Assemble Tangent Stiffness Matrix K
                // If NLGEOM is active, we use the current displacement to calculate the tangent.
                // Otherwise we always use the initial state (Linear).
                let (mut triplet, max_diag_val) = Assembler::assemble(
                    project,
                    true,
                    if self.nlgeom { Some(&u_curr) } else { None },
                )?;

                // 6. Apply Boundary Conditions (Penalty Method) to Tangent Matrix
                // For NR, we solve K * du = R.
                // Penalty on K: K_ii += Penalty
                // Penalty on R: R_i += Penalty * (u_target - u_curr_i)
                if max_diag_val > 0.0 {
                    let penalty = max_diag_val * 1.0e6;
                    for bc in &self.bcs {
                        if let Some(idx) = project.mesh.get_index_for_node_id(bc.node_id()) {
                            let global_dof = idx * 3 + bc.dof();
                            triplet.add_triplet(global_dof, global_dof, penalty);

                            let u_target = bc.value(timer, step_id);
                            residual[global_dof] += penalty * (u_target - u_curr[global_dof]);
                        }
                    }
                }

                // 7. Solve for Displacement Increment du
                let k_global: CsMat<f64> = triplet.to_csr();
                let preconditioner = DiagonalPreconditioner::new(&k_global);
                let solver: Box<dyn Solver> = match self.solver {
                    SolverType::Default | SolverType::Direct => Box::new(DirectSolver),
                    SolverType::Iterative => Box::new(IterativeSolver),
                };
                let du = solver.solve(&k_global, &residual, Some(&preconditioner), 1e-8, 1000)?;

                // 8. Update Displacement
                for i in 0..num_dofs {
                    u_curr[i] += du[i];
                }
            }

            if converged {
                current_time += increment_time;
                u_step_start.clone_from(&u_curr); // Update start of next increment

                // Create and write results for this increment
                let mut inc_solution_state = solution_state.clone();
                inc_solution_state.displacements.clone_from(&u_curr);

                let mut inc_res = IncResult::new(step_id, n_inc, "Static Step", current_time);
                let mut nodal_displacement = NodalResult::new("U", FieldType::Displacement);
                let (nodal_stress, nodal_strain) =
                    Self::calculate_stress_strain(project, &inc_solution_state)?;

                for (matrix_idx, &node_id) in project.mesh.index_to_node_id.iter().enumerate() {
                    let idx = matrix_idx * 3;
                    if idx + 2 < inc_solution_state.displacements.len() {
                        let dx = inc_solution_state.displacements[idx];
                        let dy = inc_solution_state.displacements[idx + 1];
                        let dz = inc_solution_state.displacements[idx + 2];
                        nodal_displacement.insert(node_id, vec![dx, dy, dz]);
                    }
                }
                inc_res.nodal_results.push(nodal_displacement);
                inc_res.nodal_results.push(nodal_stress);
                inc_res.nodal_results.push(nodal_strain);

                writer
                    .write_increment(&inc_res, timer)
                    .map_err(|e| FerrixError::GenericIo(std::io::Error::other(e.to_string())))?;
            } else {
                println!("Increment failed to converge. Retrying with smaller increment.");
                increment_time /= 2.;
                if increment_time < self.increment_data.min_time_increment {
                    return Err(FerrixError::ConvergenceError(
                        "Minimum time increment reached. Convergence failed.".into(),
                    ));
                }
            }
        }

        // Update the global solution state
        solution_state.displacements.clone_from(&u_step_start);

        Ok(())
    }

    /// Computes stress and strain for all nodes by averaging element contributions.
    #[allow(clippy::cast_precision_loss)]
    fn calculate_stress_strain(
        project: &Project,
        solution_state: &SolutionState,
    ) -> Result<(NodalResult, NodalResult)> {
        let mut nodal_stress = NodalResult::new("S", FieldType::Stress);
        let mut nodal_strain = NodalResult::new("E", FieldType::Strain);
        let mut node_element_count: HashMap<NodeId, usize> = HashMap::new();

        for element in project.mesh.elements.values() {
            let node_ids = element.get_node_ids();
            let mut u_el = Vec::new();
            for &node_id in node_ids {
                if let Some(idx) = project.mesh.get_index_for_node_id(node_id) {
                    let dof_start = idx * 3;
                    u_el.extend_from_slice(&solution_state.displacements[dof_start..dof_start + 3]);
                }
            }

            let material_idx = project
                .element_materials
                .get(&element.get_id())
                .ok_or_else(|| {
                    FerrixError::InvalidModelState(format!(
                        "Element {} has no material",
                        element.get_id()
                    ))
                })?;
            let material = &project.materials[*material_idx];
            let d_matrix = material.build_elastic_d_matrix(0.0)?;

            let mut avg_stress = [0.0; 6];
            let mut avg_strain = [0.0; 6];
            let integration_points = element.integration_points();
            let num_ips = integration_points.len();

            for ip in integration_points {
                let (strain, stress) =
                    element.calculate_stress_strain_at_ip(&d_matrix, &u_el, &project.mesh, &ip)?;
                for i in 0..6 {
                    avg_stress[i] += stress[i];
                    avg_strain[i] += strain[i];
                }
            }

            for i in 0..6 {
                avg_stress[i] /= num_ips as f64;
                avg_strain[i] /= num_ips as f64;
            }

            for &node_id in node_ids {
                let count = node_element_count.entry(node_id).or_insert(0);
                let current_stress = nodal_stress.data.entry(node_id).or_insert(vec![0.0; 6]);
                let current_strain = nodal_strain.data.entry(node_id).or_insert(vec![0.0; 6]);
                for i in 0..6 {
                    current_stress[i] += avg_stress[i];
                    current_strain[i] += avg_strain[i];
                }
                *count += 1;
            }
        }

        for (node_id, count) in node_element_count {
            if let Some(stress) = nodal_stress.data.get_mut(&node_id) {
                for val in stress.iter_mut() {
                    *val /= count as f64;
                }
            }
            if let Some(strain) = nodal_strain.data.get_mut(&node_id) {
                for val in strain.iter_mut() {
                    *val /= count as f64;
                }
            }
        }

        Ok((nodal_stress, nodal_strain))
    }
}
