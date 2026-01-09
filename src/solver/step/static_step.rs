use crate::solver::{
    assembler::Assembler,
    preconditioner::DiagonalPreconditioner,
    project::Project,
    solvers::{iterative::IterativeSolver, Solver},
    state::SolutionState,
    step::boundary_conds::{BoundaryCondition, Load},
};
use sprs::CsMat;
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
        loads: &[Load],
        bcs: &[BoundaryCondition],
        solution_state: &mut SolutionState,
    ) -> Result<(), Box<dyn Error>> {
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
            "System assembled. K: {}x{}, NNZ: {}. Solving...",
            k_global.rows(),
            k_global.cols(),
            k_global.nnz()
        );

        let preconditioner = DiagonalPreconditioner::new(&k_global);
        let solver = IterativeSolver;
        let delta_u = solver.solve(&k_global, &f_global, Some(&preconditioner), 1e-8, 10000)?;

        let u_norm: f64 = delta_u.iter().map(|x| x * x).sum::<f64>().sqrt();
        println!("Solution converged. Displacement Norm: {u_norm:.4e}");

        // 6. Update global solution state
        for (i, displacement) in solution_state.displacements.iter_mut().enumerate() {
            *displacement += delta_u[i];
        }

        Ok(())
    }
}