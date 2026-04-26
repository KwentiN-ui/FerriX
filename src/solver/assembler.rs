use crate::solver::error::{FerrixError, Result};
use crate::solver::project::Project;
use sprs::TriMat;

pub struct Assembler;

impl Assembler {
    /// Assembles the global stiffness matrix.
    ///
    /// # Errors
    /// Returns an error if the element stiffness matrix cannot be computed.
    pub fn assemble(project: &Project, is_symmetric: bool) -> Result<(TriMat<f64>, f64)> {
        let num_nodes = project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err(FerrixError::InvalidModelState(
                "Mesh empty or mappings not initialized".into(),
            ));
        }

        let num_dofs = num_nodes * 3;
        let mut triplet = TriMat::new((num_dofs, num_dofs));
        let mut max_diag_val: f64 = 0.0;

        for element in project.mesh.elements.values() {
            let material_index = project
                .element_materials
                .get(&element.get_id())
                .ok_or_else(|| {
                    FerrixError::InvalidModelState(format!(
                        "Element {} has no material assigned.",
                        element.get_id()
                    ))
                })?;
            let material = &project.materials[*material_index];

            let d_matrix = material.build_elastic_d_matrix();

            // Call the new high-performance method from Element
            let k_el = element.compute_stiffness(project, &d_matrix, is_symmetric)?;

            let node_ids = element.get_node_ids();
            let num_nodes_el = node_ids.len();

            for (i, node_id_i) in node_ids.iter().enumerate().take(num_nodes_el) {
                let global_index_i = project
                    .mesh
                    .get_index_for_node_id(*node_id_i)
                    .ok_or(FerrixError::NodeNotFound(*node_id_i))?;

                let start_j = if is_symmetric { i } else { 0 };

                for (j, node_id_j) in node_ids.iter().enumerate().take(num_nodes_el).skip(start_j) {
                    let global_index_j = project
                        .mesh
                        .get_index_for_node_id(*node_id_j)
                        .ok_or(FerrixError::NodeNotFound(*node_id_j))?;

                    for dof_i in 0..3 {
                        let start_dof_j = if is_symmetric && i == j { dof_i } else { 0 };

                        for dof_j in start_dof_j..3 {
                            let val = k_el[(i * 3 + dof_i, j * 3 + dof_j)];

                            if val.abs() > 1e-12 {
                                let row = global_index_i * 3 + dof_i;
                                let col = global_index_j * 3 + dof_j;

                                triplet.add_triplet(row, col, val);

                                if row == col {
                                    max_diag_val = max_diag_val.max(val.abs());
                                } else if is_symmetric {
                                    triplet.add_triplet(col, row, val);
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok((triplet, max_diag_val))
    }
}
