use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::project::Project;
use ndarray::Array2;
use sprs::TriMat;

/// Responsible for assembling the global stiffness matrix (`K`).
///
/// This struct iterates through all finite elements in the mesh, calculates
/// each element's local stiffness matrix, and assembles them into a single,
/// sparse global stiffness matrix in triplet format.
pub struct Assembler;

impl Assembler {
    /// Assembles the global stiffness matrix for the entire project.
    ///
    /// This function performs the core assembly loop:
    /// 1. Iterates over each element in the mesh.
    /// 2. Fetches the material properties for the element.
    /// 3. Calculates the element's local stiffness matrix (`k_el`).
    /// 4. Maps the local degrees of freedom to the global system.
    /// 5. Adds the `k_el` values into the global `Triplet` matrix.
    ///
    /// # Returns
    ///
    /// A `Result` containing:
    /// * A `TriMat<f64>` representing the global stiffness matrix in triplet format.
    /// * The maximum absolute value found on the diagonal, used for the penalty method.
    pub fn assemble(project: &Project) -> Result<(TriMat<f64>, f64), String> {
        let num_nodes = project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }
        let num_dofs = num_nodes * 3;
        let mut triplet = TriMat::new((num_dofs, num_dofs));
        let mut max_diag_val: f64 = 0.0;

        for element in project.mesh.elements.values() {
            let material_index = project
                .element_materials
                .get(&element.get_id())
                .ok_or(format!(
                    "Element {} has no material assigned.",
                    element.get_id()
                ))?;
            let material = &project.materials[*material_index];
            let d_matrix = material.build_elastic_d_matrix();

            let k_el = compute_element_stiffness(project, &d_matrix, element)?;
            let node_ids = element.get_node_ids();

            for (local_node_i, &global_id_i) in node_ids.iter().enumerate() {
                let global_index_i = project
                    .mesh
                    .get_index_for_node_id(global_id_i)
                    .ok_or(format!("Node {global_id_i} not found"))?;

                for (local_node_j, &global_id_j) in node_ids.iter().enumerate() {
                    let global_index_j = project
                        .mesh
                        .get_index_for_node_id(global_id_j)
                        .ok_or(format!("Node {global_id_j} not found"))?;

                    for dof_i in 0..3 {
                        for dof_j in 0..3 {
                            let val = k_el[[local_node_i * 3 + dof_i, local_node_j * 3 + dof_j]];

                            if val.abs() > 1e-12 {
                                let row = global_index_i * 3 + dof_i;
                                let col = global_index_j * 3 + dof_j;
                                triplet.add_triplet(row, col, val);

                                if row == col {
                                    max_diag_val = max_diag_val.max(val.abs());
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

fn compute_element_stiffness(
    project: &Project,
    d_mat: &Array2<f64>,
    element: &Element,
) -> Result<Array2<f64>, String> {
    let node_ids = element.get_node_ids();
    let num_nodes = node_ids.len();
    let num_dofs = num_nodes * 3;

    let mut node_coords = Array2::<f64>::zeros((3, num_nodes));
    for (i, &node_id) in node_ids.iter().enumerate() {
        let coords = project
            .mesh
            .nodes
            .get(&node_id)
            .ok_or(format!("Node {node_id} not found"))?;
        node_coords[[0, i]] = coords.x;
        node_coords[[1, i]] = coords.y;
        node_coords[[2, i]] = coords.z;
    }

    let mut k_el = Array2::<f64>::zeros((num_dofs, num_dofs));

    for gp in element.integration_points() {
        let (_, dn_local) = element.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
        let jacobian = dn_local.dot(&node_coords.t());
        let (det_j, inv_j) = invert_jacobian_3x3(&jacobian)
            .map_err(|()| format!("Singular element found with nodes: {node_ids:?}"))?;
        let dn_global = inv_j.dot(&dn_local);
        let b_mat = build_b_matrix(&dn_global, num_nodes);
        let db = d_mat.dot(&b_mat);
        let btdb = b_mat.t().dot(&db);
        k_el.scaled_add(det_j * gp.weight, &btdb);
    }

    Ok(k_el)
}

fn build_b_matrix(dn_global: &Array2<f64>, num_nodes: usize) -> Array2<f64> {
    let num_dofs = num_nodes * 3;
    let mut b = Array2::<f64>::zeros((6, num_dofs));

    for i in 0..num_nodes {
        let col_idx = i * 3;
        let d_dx = dn_global[[0, i]];
        let d_dy = dn_global[[1, i]];
        let d_dz = dn_global[[2, i]];

        b[[0, col_idx]] = d_dx;
        b[[1, col_idx + 1]] = d_dy;
        b[[2, col_idx + 2]] = d_dz;
        b[[3, col_idx]] = d_dy;
        b[[3, col_idx + 1]] = d_dx;
        b[[4, col_idx + 1]] = d_dz;
        b[[4, col_idx + 2]] = d_dy;
        b[[5, col_idx]] = d_dz;
        b[[5, col_idx + 2]] = d_dx;
    }
    b
}

fn invert_jacobian_3x3(m: &Array2<f64>) -> Result<(f64, Array2<f64>), ()> {
    let det = m[[0, 0]] * (m[[1, 1]] * m[[2, 2]] - m[[2, 1]] * m[[1, 2]])
        - m[[0, 1]] * (m[[1, 0]] * m[[2, 2]] - m[[1, 2]] * m[[2, 0]])
        + m[[0, 2]] * (m[[1, 0]] * m[[2, 1]] - m[[1, 1]] * m[[2, 0]]);

    if det.abs() < 1e-14 {
        return Err(());
    }

    let inv_det = 1.0 / det;
    let mut inv = Array2::<f64>::zeros((3, 3));

    inv[[0, 0]] = (m[[1, 1]] * m[[2, 2]] - m[[2, 1]] * m[[1, 2]]) * inv_det;
    inv[[0, 1]] = (m[[0, 2]] * m[[2, 1]] - m[[0, 1]] * m[[2, 2]]) * inv_det;
    inv[[0, 2]] = (m[[0, 1]] * m[[1, 2]] - m[[0, 2]] * m[[1, 1]]) * inv_det;

    inv[[1, 0]] = (m[[1, 2]] * m[[2, 0]] - m[[1, 0]] * m[[2, 2]]) * inv_det;
    inv[[1, 1]] = (m[[0, 0]] * m[[2, 2]] - m[[0, 2]] * m[[2, 0]]) * inv_det;
    inv[[1, 2]] = (m[[1, 0]] * m[[0, 2]] - m[[0, 0]] * m[[1, 2]]) * inv_det;

    inv[[2, 0]] = (m[[1, 0]] * m[[2, 1]] - m[[2, 0]] * m[[1, 1]]) * inv_det;
    inv[[2, 1]] = (m[[2, 0]] * m[[0, 1]] - m[[0, 0]] * m[[2, 1]]) * inv_det;
    inv[[2, 2]] = (m[[0, 0]] * m[[1, 1]] - m[[1, 0]] * m[[0, 1]]) * inv_det;

    Ok((det, inv))
}
