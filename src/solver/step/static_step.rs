use crate::solver::assembler::Assembler;
use crate::solver::ids::NodeId;
use crate::solver::increment::IncrementData;
use crate::solver::io::writer::ResultWriter;
use crate::solver::preconditioner::DiagonalPreconditioner;
use crate::solver::project::Project;
use crate::solver::results::{FieldType, IncResult, NodalResult};
use crate::solver::solvers::SolverType;
use crate::solver::solvers::{Solver, direct::DirectSolver, iterative::IterativeSolver};
use crate::solver::state::SolutionState;
use crate::solver::time::SolverTime;
use sprs::CsMat;
use std::collections::HashMap;
use std::error::Error;

#[derive(Debug, Clone)]
pub struct StaticStep {
    pub solver: SolverType,
    pub increment_data: IncrementData,
}

impl StaticStep {
    pub fn solve(
        &self,
        step_id: usize,
        project: &Project,
        solution_state: &mut SolutionState,
        writer: &dyn ResultWriter,
        timer: &mut SolverTime,
    ) -> Result<(), Box<dyn Error>> {
        println!("--- Step {step_id}: StaticStep ---");

        let mut current_time = 0.0;
        let mut n_inc = 0;
        let mut increment_time = self.increment_data.initial_time_increment;
        let mut step_displacement_increment = vec![0.0; solution_state.displacements.len()];

        while current_time < self.increment_data.time_period {
            n_inc += 1;
            if n_inc > self.increment_data.max_iterations {
                return Err(format!(
                    "Maximum increment count of {} exceeded.",
                    self.increment_data.max_iterations,
                )
                .into());
            }

            if current_time + increment_time > self.increment_data.time_period {
                increment_time = self.increment_data.time_period - current_time;
            }

            timer.new_increment(increment_time);

            println!(
                "Increment {n_inc} | Step Time: {:.4e} | Increment Size: {increment_time:.4e}",
                timer.local_time()
            );

            match self.next_increment(project, timer) {
                Ok(delta_u) => {
                    current_time += increment_time;

                    // Accumulate displacements for the step
                    for (i, val) in delta_u.iter().enumerate() {
                        step_displacement_increment[i] += val;
                    }

                    // Create a temporary solution state for this increment's results
                    let mut inc_solution_state = solution_state.clone();
                    for (i, val) in step_displacement_increment.iter().enumerate() {
                        inc_solution_state.displacements[i] += val;
                    }

                    // Create and write results for this increment
                    let mut inc_res = IncResult::new(step_id, n_inc, "Static Step", current_time);

                    let mut nodal_displacement = NodalResult::new("U", FieldType::Displacement);
                    let (nodal_stress, nodal_strain) =
                        Self::calculate_stress_strain(project, &inc_solution_state);

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

                    writer.write_increment(&inc_res)?;
                }
                Err(e) => {
                    println!("Increment failed: {e}. Retrying with smaller increment.");
                    increment_time /= 2.;
                    if increment_time < self.increment_data.min_time_increment {
                        return Err("Minimum time increment reached. Convergence failed.".into());
                    }
                }
            }
        }

        // Update the global solution state once at the end of the step
        for (i, val) in step_displacement_increment.iter().enumerate() {
            solution_state.displacements[i] += val;
        }

        Ok(())
    }

    fn next_increment(&self, project: &Project, timer: &SolverTime) -> Result<Vec<f64>, String> {
        // 1. Setup
        let num_nodes = project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }
        let num_dofs = num_nodes * 3;

        // Init Force Vector F (for the increment)
        let mut f_inc = vec![0.0; num_dofs];

        // 2. Add loads to F
        for load in &project.loads {
            if let Some(idx) = project.mesh.get_index_for_node_id(load.node_id()) {
                let global_dof = idx * 3 + load.dof();
                if global_dof < num_dofs {
                    f_inc[global_dof] += load.value(timer);
                }
            } else {
                eprintln!("Warning: Load on unknown node {}", load.node_id());
            }
        }

        // 3. Assemble stiffness matrix
        let (mut triplet, max_diag_val) = Assembler::assemble(project)?;

        // 4. Apply boundary conditions (Penalty Method)
        if max_diag_val > 0.0 {
            let penalty = max_diag_val * 1.0e6;
            for bc in &project.bcs {
                if let Some(idx) = project.mesh.get_index_for_node_id(bc.node_id()) {
                    let global_dof = idx * 3 + bc.dof();
                    if global_dof < num_dofs {
                        triplet.add_triplet(global_dof, global_dof, penalty);
                        f_inc[global_dof] += penalty * bc.value(timer);
                    }
                }
            }
        }

        // 5. Conversion & Solving for displacement increment
        let k_global: CsMat<f64> = triplet.to_csr();
        let preconditioner = DiagonalPreconditioner::new(&k_global);
        let solver: Box<dyn Solver> = match self.solver {
            crate::solver::solvers::SolverType::Default
            | crate::solver::solvers::SolverType::Direct => Box::new(DirectSolver),
            crate::solver::solvers::SolverType::Iterative => Box::new(IterativeSolver),
        };
        let delta_u = solver.solve(&k_global, &f_inc, Some(&preconditioner), 1e-8, 10000)?;
        Ok(delta_u)
    }

    #[allow(clippy::cast_precision_loss)]
    fn calculate_stress_strain(
        project: &Project,
        solution_state: &SolutionState,
    ) -> (NodalResult, NodalResult) {
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

            let material = &project.materials[project.element_materials[&element.get_id()]];
            let d_matrix = material.build_elastic_d_matrix();

            let mut avg_stress = [0.0; 6];
            let mut avg_strain = [0.0; 6];
            let integration_points = element.integration_points();
            let num_ips = integration_points.len();

            for ip in integration_points {
                let (strain, stress) =
                    element.calculate_stress_strain_at_ip(&d_matrix, &u_el, &project.mesh, &ip);
                for i in 0..6 {
                    avg_stress[i] += stress[i];
                    avg_strain[i] += strain[i];
                }
            }

            for i in 0..6 {
                avg_stress[i] /= num_ips as f64;
                avg_strain[i] /= num_ips as f64;
            }

            // This is a simple averaging scheme. A more sophisticated approach would be to
            // extrapolate from Gauss points to nodes.
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

        (nodal_stress, nodal_strain)
    }
}
