use crate::solver::ids::NodeId;
use crate::solver::solvers::SolverType;
use crate::solver::solvers::direct::DirectSolver;
use crate::solver::{
    assembler::Assembler,
    preconditioner::DiagonalPreconditioner,
    project::Project,
    results::{FieldType, NodalResult, StepResult},
    solvers::{Solver, iterative::IterativeSolver},
    state::SolutionState,
    step::boundary_conds::{BoundaryCondition, Load},
};
use sprs::CsMat;
use std::collections::HashMap;
use std::error::Error;

/// Represents a static analysis step in the FEA simulation.
///
/// This struct holds the state and logic required to perform a linear static analysis,
/// which solves the equation `K * u = F`, where `K` is the global stiffness matrix,
/// `u` is the displacement vector, and `F` is the external force vector.
#[derive(Debug, Clone)]
pub struct StaticStep {
    project: Box<Project>,
}

impl StaticStep {
    pub fn new(project: Box<Project>) -> Self {
        Self { project }
    }

    /// Executes the static analysis step.
    ///
    /// This is the main entry point for running the simulation for this step. It orchestrates
    /// the entire process, including:
    /// 1. Assembling the global stiffness matrix `K`.
    /// 2. Constructing the global force vector `F`.
    /// 3. Applying boundary conditions to the system of equations.
    /// 4. Solving for the incremental displacement vector `delta_u`.
    /// 5. Updating the global `SolutionState` with the incremental displacements.
    ///
    /// # Arguments
    ///
    /// * `loads` - A slice of `Load` structs representing the external forces for this increment.
    /// * `bcs` - A slice of `BoundaryCondition` structs representing the constraints for this increment.
    /// * `solution_state` - A mutable reference to the global `SolutionState`.
    pub fn compute(
        &mut self,
        step_id: usize,
        loads: &[Load],
        bcs: &[BoundaryCondition],
        solution_state: &mut SolutionState,
        solver: &SolverType,
    ) -> Result<StepResult, Box<dyn Error>> {
        // 1. Setup
        println!("Constructing global stiffness matrix");
        let num_nodes = self.project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }
        let num_dofs = num_nodes * 3;

        // Init Force Vector F
        let mut f_global = vec![0.0; num_dofs];

        // 2. Add loads to F
        for load in loads {
            if let Some(idx) = self.project.mesh.get_index_for_node_id(load.node_id) {
                let global_dof = idx * 3 + load.dof;
                if global_dof < num_dofs {
                    f_global[global_dof] += load.value;
                }
            } else {
                eprintln!("Warning: Load on unknown node {}", load.node_id);
            }
        }

        // 3. Assemble stiffness matrix
        let (mut triplet, max_diag_val) = Assembler::assemble(&self.project)?;

        // 4. Apply boundary conditions (Penalty Method)
        if max_diag_val > 0.0 {
            let penalty = max_diag_val * 1.0e6;
            for bc in bcs {
                if let Some(idx) = self.project.mesh.get_index_for_node_id(bc.node_id) {
                    let global_dof = idx * 3 + bc.dof;
                    if global_dof < num_dofs {
                        triplet.add_triplet(global_dof, global_dof, penalty);
                        f_global[global_dof] += penalty * bc.value;
                    }
                }
            }
        }

        // 5. Conversion & Solving
        let k_global: CsMat<f64> = triplet.to_csr();

        println!(
            "System assembled. K: {}x{}, NNZ: {}",
            k_global.rows(),
            k_global.cols(),
            k_global.nnz()
        );

        let preconditioner = DiagonalPreconditioner::new(&k_global);
        // let solver = IterativeSolver;
        let solver: Box<dyn Solver> = match solver {
            SolverType::Default | SolverType::Direct => Box::new(DirectSolver),
            SolverType::Iterative => Box::new(IterativeSolver),
        };
        let delta_u = solver.solve(&k_global, &f_global, Some(&preconditioner), 1e-8, 10000)?;

        let u_norm: f64 = delta_u.iter().map(|x| x * x).sum::<f64>().sqrt();
        println!("Solution found. Displacement Norm: {u_norm:.4e}");

        // 6. Update global solution state
        for (i, displacement) in solution_state.displacements.iter_mut().enumerate() {
            *displacement += delta_u[i];
        }

        // --- Create results for this step ---
        let mut step_res = StepResult::new(step_id, "Static Step", 1.0);

        // Nodal results
        let mut nodal_displacement = NodalResult::new("U", FieldType::Displacement);
        let (nodal_stress, nodal_strain) = self.calculate_stress_strain(solution_state);

        for (matrix_idx, &node_id) in self.project.mesh.index_to_node_id.iter().enumerate() {
            let idx = matrix_idx * 3;
            if idx + 2 < solution_state.displacements.len() {
                let dx = solution_state.displacements[idx];
                let dy = solution_state.displacements[idx + 1];
                let dz = solution_state.displacements[idx + 2];
                nodal_displacement.insert(node_id, vec![dx, dy, dz]);
            }
        }
        step_res.nodal_results.push(nodal_displacement);
        step_res.nodal_results.push(nodal_stress);
        step_res.nodal_results.push(nodal_strain);

        Ok(step_res)
    }

    #[allow(clippy::cast_precision_loss)]
    fn calculate_stress_strain(
        &self,
        solution_state: &SolutionState,
    ) -> (NodalResult, NodalResult) {
        let mut nodal_stress = NodalResult::new("S", FieldType::Stress);
        let mut nodal_strain = NodalResult::new("E", FieldType::Strain);
        let mut node_element_count: HashMap<NodeId, usize> = HashMap::new();

        for element in self.project.mesh.elements.values() {
            let node_ids = element.get_node_ids();
            let mut u_el = Vec::new();
            for &node_id in node_ids {
                if let Some(idx) = self.project.mesh.get_index_for_node_id(node_id) {
                    let dof_start = idx * 3;
                    u_el.extend_from_slice(&solution_state.displacements[dof_start..dof_start + 3]);
                }
            }

            let material =
                &self.project.materials[self.project.element_materials[&element.get_id()]];
            let d_matrix = material.build_elastic_d_matrix();

            let mut avg_stress = [0.0; 6];
            let mut avg_strain = [0.0; 6];
            let integration_points = element.integration_points();
            let num_ips = integration_points.len();

            for ip in integration_points {
                let (strain, stress) = element.calculate_stress_strain_at_ip(
                    &d_matrix,
                    &u_el,
                    &self.project.mesh,
                    &ip,
                );
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
