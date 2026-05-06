//! Global system assembly logic.
//!
//! This module provides the `Assembler`, which iterates over the mesh elements,
//! computes their local stiffness and internal force contributions, and
//! assembles them into the global sparse system.

use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::ElementId;
use crate::solver::material::ElementMaterialState;
use crate::solver::project::Project;
use rayon::prelude::*;
use sprs::TriMat;
use std::collections::HashMap;

/// Type definition for material states mapping elements to their integration point variables.
pub type MaterialStates = HashMap<ElementId, ElementMaterialState>;

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
    /// * `u_global` - Optional current displacement field (used for non-linear stiffness).
    /// * `t_initial` - Optional nodal temperatures at start of analysis.
    /// * `t_current` - Optional nodal temperatures at current time.
    /// * `material_states_old` - Optional SDVs at the start of the increment.
    /// * `dtime` - Time increment.
    ///
    /// # Errors
    /// Returns `FerrixError` if an element has no material or if node mappings are missing.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble(
        project: &Project,
        u_global: Option<&[f64]>,
        t_initial: Option<&[f64]>,
        t_current: Option<&[f64]>,
        material_states_old: Option<&MaterialStates>,
        dtime: f64,
    ) -> Result<(TriMat<f64>, f64, MaterialStates)> {
        let num_nodes = project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err(FerrixError::InvalidModelState(
                "Mesh empty or mappings not initialized".into(),
            ));
        }

        let num_dofs = num_nodes * 3;
        let mut material_states_new = HashMap::new();

        let (triplets, global_max_diag, updated_states): (Vec<_>, f64, Vec<_>) = project
            .mesh
            .elements
            .values()
            .par_bridge()
            .map(|element| {
                let elem_id = element.get_id();
                let material_index = project.element_materials.get(&elem_id).ok_or_else(|| {
                    FerrixError::InvalidModelState(format!(
                        "Element {elem_id} has no material assigned."
                    ))
                })?;
                let material = &project.materials[*material_index];

                let node_ids = element.get_node_ids();

                let extract_temps = |t_glob: Option<&[f64]>| -> Result<Option<Vec<f64>>> {
                    if let Some(t) = t_glob {
                        let mut t_vec = Vec::with_capacity(node_ids.len());
                        for &node_id in node_ids {
                            let global_idx = project
                                .mesh
                                .get_index_for_node_id(node_id)
                                .ok_or(FerrixError::NodeNotFound(node_id))?;
                            t_vec.push(t[global_idx]);
                        }
                        Ok(Some(t_vec))
                    } else {
                        Ok(None)
                    }
                };

                let temps_init = extract_temps(t_initial)?;
                let temps_curr = extract_temps(t_current)?;

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

                let elem_states_old = material_states_old.and_then(|m| m.get(&elem_id));

                let (k_el, updated_states) = element.compute_stiffness_sdv(
                    project,
                    material.as_ref(),
                    u_el.as_deref(),
                    temps_init.as_deref(),
                    temps_curr.as_deref(),
                    elem_states_old,
                    dtime,
                )?;

                let mut local_triplets = Vec::new();
                let mut local_max_diag: f64 = 0.0;
                let num_nodes_el = node_ids.len();

                for (i, node_id_i) in node_ids.iter().enumerate().take(num_nodes_el) {
                    let global_index_i = project
                        .mesh
                        .get_index_for_node_id(*node_id_i)
                        .ok_or(FerrixError::NodeNotFound(*node_id_i))?;

                    for (j, node_id_j) in node_ids.iter().enumerate().take(num_nodes_el) {
                        let global_index_j = project
                            .mesh
                            .get_index_for_node_id(*node_id_j)
                            .ok_or(FerrixError::NodeNotFound(*node_id_j))?;

                        for dof_i in 0..3 {
                            for dof_j in 0..3 {
                                let val = k_el[(i * 3 + dof_i, j * 3 + dof_j)];

                                if val.abs() > 1e-12 {
                                    let row = global_index_i * 3 + dof_i;
                                    let col = global_index_j * 3 + dof_j;

                                    local_triplets.push((row, col, val));

                                    if row == col {
                                        local_max_diag = local_max_diag.max(val.abs());
                                    }
                                }
                            }
                        }
                    }
                }

                Ok((local_triplets, local_max_diag, (elem_id, updated_states)))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(
                (Vec::new(), 0.0, Vec::new()),
                |(mut acc_t, mut acc_m, mut acc_s), (t, m, s)| {
                    acc_t.extend(t);
                    acc_m = acc_m.max(m);
                    acc_s.push(s);
                    (acc_t, acc_m, acc_s)
                },
            );

        let mut triplet = TriMat::new((num_dofs, num_dofs));
        for (row, col, val) in triplets {
            triplet.add_triplet(row, col, val);
        }

        for (elem_id, state) in updated_states {
            if let Some(s) = state {
                material_states_new.insert(elem_id, s);
            }
        }

        Ok((triplet, global_max_diag, material_states_new))
    }

    /// Assembles the global internal force vector.
    ///
    /// # Errors
    /// Returns an error if the element internal forces cannot be computed.
    pub fn assemble_internal_force(
        project: &Project,
        u_global: &[f64],
        t_initial: &[f64],
        t_current: &[f64],
        material_states_old: Option<&MaterialStates>,
        dtime: f64,
        u_conf: Option<&[f64]>,
    ) -> Result<(Vec<f64>, MaterialStates)> {
        let num_nodes = project.mesh.nodes.len();
        let num_dofs = num_nodes * 3;
        let mut f_int_global = vec![0.0; num_dofs];
        let mut material_states_new = HashMap::new();

        let (f_ints, updated_states): (Vec<_>, Vec<_>) = project
            .mesh
            .elements
            .values()
            .par_bridge()
            .map(|element| {
                let elem_id = element.get_id();
                let material_index = project.element_materials.get(&elem_id).ok_or_else(|| {
                    FerrixError::InvalidModelState(format!(
                        "Element {elem_id} has no material assigned."
                    ))
                })?;
                let material = &project.materials[*material_index];

                let node_ids = element.get_node_ids();
                let mut node_temps_init = Vec::with_capacity(node_ids.len());
                let mut node_temps_curr = Vec::with_capacity(node_ids.len());
                for &node_id in node_ids {
                    let global_idx = project
                        .mesh
                        .get_index_for_node_id(node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    node_temps_init.push(t_initial[global_idx]);
                    node_temps_curr.push(t_current[global_idx]);
                }

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

                let elem_states_old = material_states_old.and_then(|m| m.get(&elem_id));

                let (f_int_el, updated_states) = element.compute_internal_force_sdv(
                    &project.mesh,
                    material.as_ref(),
                    &u_el,
                    &node_temps_init,
                    &node_temps_curr,
                    elem_states_old,
                    dtime,
                    u_conf_el.as_deref(),
                )?;

                let mut local_f_int = Vec::new();
                for (i, &node_id) in node_ids.iter().enumerate() {
                    let global_idx = project
                        .mesh
                        .get_index_for_node_id(node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    for dof in 0..3 {
                        local_f_int.push((global_idx * 3 + dof, f_int_el[i * 3 + dof]));
                    }
                }

                Ok((local_f_int, (elem_id, updated_states)))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .unzip();

        for local_f_int in f_ints {
            for (idx, val) in local_f_int {
                f_int_global[idx] += val;
            }
        }

        for (elem_id, state) in updated_states {
            if let Some(s) = state {
                material_states_new.insert(elem_id, s);
            }
        }

        Ok((f_int_global, material_states_new))
    }

    /// Assembles the global thermal force vector.
    ///
    /// # Errors
    /// Returns an error if any element thermal force cannot be computed.
    pub fn assemble_thermal_force(
        project: &Project,
        t_initial: &[f64],
        t_current: &[f64],
    ) -> Result<Vec<f64>> {
        let num_nodes = project.mesh.nodes.len();
        let num_dofs = num_nodes * 3;
        let mut f_th_global = vec![0.0; num_dofs];

        let f_ths: Vec<_> = project
            .mesh
            .elements
            .values()
            .par_bridge()
            .map(|element| {
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

                let node_ids = element.get_node_ids();
                let mut node_temps_init = Vec::with_capacity(node_ids.len());
                let mut node_temps_curr = Vec::with_capacity(node_ids.len());
                for &node_id in node_ids {
                    let global_idx = project
                        .mesh
                        .get_index_for_node_id(node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    node_temps_init.push(t_initial[global_idx]);
                    node_temps_curr.push(t_current[global_idx]);
                }

                let f_th_el = element.compute_thermal_force(
                    project,
                    material.as_ref(),
                    &node_temps_init,
                    &node_temps_curr,
                )?;

                let mut local_f_th = Vec::new();
                for (i, &node_id) in node_ids.iter().enumerate() {
                    let global_idx = project
                        .mesh
                        .get_index_for_node_id(node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    for dof in 0..3 {
                        local_f_th.push((global_idx * 3 + dof, f_th_el[i * 3 + dof]));
                    }
                }

                Ok(local_f_th)
            })
            .collect::<Result<Vec<_>>>()?;

        for local_f_th in f_ths {
            for (idx, val) in local_f_th {
                f_th_global[idx] += val;
            }
        }

        Ok(f_th_global)
    }
}
