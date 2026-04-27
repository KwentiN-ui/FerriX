//! Global system assembly logic.
//!
//! This module provides the `Assembler`, which iterates over the mesh elements,
//! computes their local stiffness and internal force contributions, and
//! assembles them into the global sparse system.

use crate::solver::error::{FerrixError, Result};
use crate::solver::project::Project;
use sprs::TriMat;

/// A utility for assembling global matrices and vectors from element contributions.
pub struct Assembler;

impl Assembler {
    /// Assembles the global stiffness matrix (K-matrix).
    ///
    /// Iterates through all elements in the project, computes their local stiffness
    /// matrices based on material properties, and maps them to the global system of equations.
    ///
    /// # Arguments
    /// * `project` - The FEA project containing mesh and materials.
    /// * `is_symmetric` - Whether to assume and enforce matrix symmetry.
    /// * `u_global` - Optional current displacement field (used for non-linear stiffness).
    ///
    /// # Errors
    /// Returns `FerrixError` if an element has no material or if node mappings are missing.
    pub fn assemble(
        project: &Project,
        is_symmetric: bool,
        u_global: Option<&[f64]>,
    ) -> Result<(TriMat<f64>, f64)> {
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

            let d_matrix = material.build_elastic_d_matrix(0.0)?;

            let node_ids = element.get_node_ids();

            // Extract u_el if u_global is provided
            let u_el = if let Some(u_glob) = u_global {
                let mut u_vec = Vec::with_capacity(node_ids.len() * 3);
                for &node_id in node_ids {
                    let global_idx = project
                        .mesh
                        .get_index_for_node_id(node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    u_vec.extend_from_slice(&u_glob[global_idx * 3..global_idx * 3 + 3]);
                }
                Some(u_vec)
            } else {
                None
            };

            // Call compute_stiffness with u_el
            let k_el =
                element.compute_stiffness(project, &d_matrix, is_symmetric, u_el.as_deref())?;

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

    /// Assembles the global internal force vector.
    ///
    /// # Errors
    /// Returns an error if the element internal forces cannot be computed.
    pub fn assemble_internal_force(
        project: &Project,
        u_global: &[f64],
        u_conf: Option<&[f64]>,
    ) -> Result<Vec<f64>> {
        let num_nodes = project.mesh.nodes.len();
        let num_dofs = num_nodes * 3;
        let mut f_int_global = vec![0.0; num_dofs];

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
            let d_matrix = material.build_elastic_d_matrix(0.0)?;

            let node_ids = element.get_node_ids();
            let mut u_el = Vec::with_capacity(node_ids.len() * 3);
            let mut u_conf_el = if u_conf.is_some() {
                Some(Vec::with_capacity(node_ids.len() * 3))
            } else {
                None
            };

            for &node_id in node_ids {
                let global_idx = project
                    .mesh
                    .get_index_for_node_id(node_id)
                    .ok_or(FerrixError::NodeNotFound(node_id))?;
                u_el.extend_from_slice(&u_global[global_idx * 3..global_idx * 3 + 3]);
                if let Some(uc_el) = &mut u_conf_el {
                    if let Some(uc_glob) = u_conf {
                        uc_el.extend_from_slice(&uc_glob[global_idx * 3..global_idx * 3 + 3]);
                    }
                }
            }

            let f_int_el = element.compute_internal_force(
                &project.mesh,
                &u_el,
                &d_matrix,
                u_conf_el.as_deref(),
            )?;

            for (i, _) in node_ids.iter().enumerate() {
                let global_idx = project
                    .mesh
                    .get_index_for_node_id(node_ids[i])
                    .ok_or(FerrixError::NodeNotFound(node_ids[i]))?;
                for dof in 0..3 {
                    f_int_global[global_idx * 3 + dof] += f_int_el[i * 3 + dof];
                }
            }
        }

        Ok(f_int_global)
    }
}
