use crate::solver::assembler::Assembler;
#[allow(unused_imports)]
use crate::solver::ids::{BoundaryConditionId, LoadId, NodeId};
use crate::solver::{
    project::Project,
    results::{FieldType, NodalResult, StepResult},
    solver::{IterativeSolver, Solver},
    step::boundary_conds::{BoundaryCondition, Load},
};
use sprs::CsMat;
use std::error::Error;

#[derive(Debug, Clone)]
pub struct StaticStep {
    project: Box<Project>,
}

impl StaticStep {
    pub fn new(project: Box<Project>) -> Self {
        Self { project }
    }

    pub fn compute(
        &mut self,
        step_id: usize,
        loads: &[Load],
        bcs: &[BoundaryCondition],
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
            "System assembled. K: {}x{}, NNZ: {}. Solving...",
            k_global.rows(),
            k_global.cols(),
            k_global.nnz()
        );

        let solver = IterativeSolver;
        let u = solver.solve(&k_global, &f_global, 1e-8, 10000)?;

        let u_norm: f64 = u.iter().map(|x| x * x).sum::<f64>().sqrt();
        println!("Solution converged. Displacement Norm: {u_norm:.4e}");

        let mut displacement_field = NodalResult::new("U", FieldType::Displacement);
        for (matrix_idx, &node_id) in self.project.mesh.index_to_node_id.iter().enumerate() {
            let idx = matrix_idx * 3;
            if idx + 2 < u.len() {
                let dx = u[idx];
                let dy = u[idx + 1];
                let dz = u[idx + 2];
                displacement_field.insert(node_id, vec![dx, dy, dz]);
            }
        }

        let mut step_res = StepResult::new(step_id, "Static Step", 1.);
        step_res.nodal_results.push(displacement_field);

        Ok(step_res)
    }
}