//! Static stress analysis step.
//!
//! This module implements the solver for static FEA problems, supporting both
//! linear and non-linear (NLGEOM) analyses using the Newton-Raphson method.

use crate::solver::assembler::Assembler;
use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::NodeId;
use crate::solver::increment::IncrementData;
use crate::solver::io::writer::ResultWriter;
use crate::solver::material::MaterialPointState;
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

        // Start of increment states
        let mut states_inc_start = solution_state.material_states.clone();

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

            // material_states_curr will store the SDVs calculated during the iteration
            let mut material_states_curr = states_inc_start.clone();

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

                // Add Thermal Forces
                let f_th = Assembler::assemble_thermal_force(
                    project,
                    &solution_state.initial_temperatures,
                    &solution_state.temperatures,
                )?;
                for i in 0..num_dofs {
                    f_ext[i] += f_th[i];
                }

                // 2. Calculate Internal Forces F_int and Updated States
                let u_conf = if self.nlgeom {
                    Some(u_curr.as_slice())
                } else {
                    None
                };

                let (f_int, updated_states) = Assembler::assemble_internal_force(
                    project,
                    &u_curr,
                    &solution_state.initial_temperatures,
                    &solution_state.temperatures,
                    Some(&states_inc_start),
                    increment_time,
                    u_conf,
                )?;
                material_states_curr = updated_states;

                // 3. Calculate Residual R = F_ext - F_int
                let mut residual = vec![0.0; num_dofs];
                for i in 0..num_dofs {
                    residual[i] = f_ext[i] - f_int[i];
                }

                // 4. Check Convergence
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

                if rel_res < 1e-3 && iter > 0 {
                    converged = true;
                    break;
                }

                // 5. Assemble Tangent Stiffness Matrix K
                let (mut triplet, max_diag_val, _) = Assembler::assemble(
                    project,
                    if self.nlgeom { Some(&u_curr) } else { None },
                    Some(&solution_state.initial_temperatures),
                    Some(&solution_state.temperatures),
                    Some(&states_inc_start),
                    increment_time,
                )?;

                // 6. Apply Boundary Conditions
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

                // 7. Solve
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
                u_step_start.clone_from(&u_curr);
                states_inc_start.clone_from(&material_states_curr); // Commit states for next increment

                // Create and write results
                let mut inc_solution_state = solution_state.clone();
                inc_solution_state.displacements.clone_from(&u_curr);
                inc_solution_state
                    .material_states
                    .clone_from(&material_states_curr);

                let mut inc_res = IncResult::new(step_id, n_inc, "Static Step", current_time);
                let mut nodal_displacement = NodalResult::new("U", FieldType::Displacement);
                let mut nodal_temperature = NodalResult::new("NT", FieldType::Temperature);
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
                    if matrix_idx < inc_solution_state.temperatures.len() {
                        nodal_temperature
                            .insert(node_id, vec![inc_solution_state.temperatures[matrix_idx]]);
                    }
                }
                inc_res.nodal_results.push(nodal_displacement);
                inc_res.nodal_results.push(nodal_temperature);
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
        solution_state.material_states.clone_from(&states_inc_start);

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
            let elem_id = element.get_id();
            let node_ids = element.get_node_ids();
            let mut u_el = Vec::new();
            let mut t_init_el = Vec::new();
            let mut t_curr_el = Vec::new();
            for &node_id in node_ids {
                if let Some(idx) = project.mesh.get_index_for_node_id(node_id) {
                    let dof_start = idx * 3;
                    u_el.extend_from_slice(&solution_state.displacements[dof_start..dof_start + 3]);
                    t_init_el.push(solution_state.initial_temperatures[idx]);
                    t_curr_el.push(solution_state.temperatures[idx]);
                }
            }

            let material_idx = project.element_materials.get(&elem_id).ok_or_else(|| {
                FerrixError::InvalidModelState(format!("Element {elem_id} has no material"))
            })?;
            let material = &project.materials[*material_idx];

            let t_avg = t_curr_el.iter().sum::<f64>() / t_curr_el.len() as f64;
            let d_matrix = material.build_elastic_d_matrix(t_avg)?;

            let node_local_coords = element.node_local_coords();
            let num_ips = element.integration_points().len();

            // Calculate centroid for fallback
            let mut centroid = [0.0; 3];
            for p in &node_local_coords {
                centroid[0] += p[0];
                centroid[1] += p[1];
                centroid[2] += p[2];
            }
            let n_nodes = node_local_coords.len() as f64;
            let centroid = [
                centroid[0] / n_nodes,
                centroid[1] / n_nodes,
                centroid[2] / n_nodes,
            ];

            for (i, &node_id) in node_ids.iter().enumerate() {
                let local_pos = node_local_coords[i];
                let res = element.calculate_stress_strain_at_local_coords(
                    &d_matrix,
                    &u_el,
                    &project.mesh,
                    local_pos[0],
                    local_pos[1],
                    local_pos[2],
                );

                let (strain_total, _stress) = match res {
                    Ok(s) => s,
                    Err(FerrixError::NumericalError(ref msg))
                        if msg.contains("Singular Jacobian") =>
                    {
                        // Fallback: shift slightly towards element center
                        let eps = 1e-3;
                        let shifted = [
                            local_pos[0] * (1.0 - eps) + centroid[0] * eps,
                            local_pos[1] * (1.0 - eps) + centroid[1] * eps,
                            local_pos[2] * (1.0 - eps) + centroid[2] * eps,
                        ];
                        element.calculate_stress_strain_at_local_coords(
                            &d_matrix,
                            &u_el,
                            &project.mesh,
                            shifted[0],
                            shifted[1],
                            shifted[2],
                        )?
                    }
                    Err(e) => return Err(e),
                };

                let mut strain_mech = strain_total;
                let t_curr_node = t_curr_el[i];
                let t_init_node = t_init_el[i];

                if let Some(alpha) = material.thermal_expansion(t_curr_node) {
                    let t_ref = material.reference_temperature();
                    let th_strain_curr = alpha * (t_curr_node - t_ref);

                    let alpha_init = material.thermal_expansion(t_init_node).unwrap_or(alpha);
                    let th_strain_init = alpha_init * (t_init_node - t_ref);

                    let delta_th_strain = th_strain_curr - th_strain_init;

                    strain_mech[0] -= delta_th_strain;
                    strain_mech[1] -= delta_th_strain;
                    strain_mech[2] -= delta_th_strain;
                }

                let mut avg_state = MaterialPointState::default();
                let num_sdvs = material.num_state_variables();
                if num_sdvs > 0 {
                    avg_state.variables = vec![0.0; num_sdvs];
                    let elem_state = solution_state.material_states.get(&elem_id);
                    if let Some(es) = elem_state {
                        for ip_idx in 0..num_ips {
                            for sdv_idx in 0..num_sdvs {
                                avg_state[sdv_idx] += es[ip_idx][sdv_idx] / num_ips as f64;
                            }
                        }
                    }
                }

                // If we used fallback, we should ideally recompute stress with updated strain_mech
                // but update_state for linear elasticity is just stress = D * strain.
                // However, material.update_state might be non-linear.
                // For now, if we had to fallback, the stress returned by calculate_stress_strain_at_local_coords
                // was using total strain. Let's re-calculate stress from strain_mech.
                let (_, stress_mech, _) =
                    material.update_state(t_curr_node, &strain_mech, &avg_state, 0.0)?;

                let current_stress = nodal_stress.data.entry(node_id).or_insert(vec![0.0; 6]);
                let current_strain = nodal_strain.data.entry(node_id).or_insert(vec![0.0; 6]);
                for j in 0..6 {
                    current_stress[j] += stress_mech[j];
                    current_strain[j] += strain_mech[j];
                }
                *node_element_count.entry(node_id).or_insert(0) += 1;
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
