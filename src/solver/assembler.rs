use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::project::Project;
use nalgebra::{DMatrix, SMatrix};
use sprs::TriMat;

pub struct Assembler;

impl Assembler {
    pub fn assemble(project: &Project) -> Result<(TriMat<f64>, f64), String> {
        let num_nodes = project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }
        // Each node has 3 Degrees of Freedom (DOFs) in 3D: u_x, u_y, u_z
        let num_dofs = num_nodes * 3;
        let mut triplet = TriMat::new((num_dofs, num_dofs));
        let mut max_diag_val: f64 = 0.0;

        for element in project.mesh.elements.values() {
            let material_index =
                project
                    .element_materials
                    .get(&element.get_id())
                    .ok_or(format!(
                        "Element {} has no material assigned.",
                        element.get_id()
                    ))?;
            let material = &project.materials[*material_index];

            // 1. Material Law (D-matrix): Defines the stress-strain relationship (Hooke's Law)
            let d_matrix = material.build_elastic_d_matrix();

            // 2. Compute local stiffness matrix k_el (the integral of B^T * D * B)
            let k_el = compute_element_stiffness(project, &d_matrix, element)?;
            let node_ids = element.get_node_ids();

            // 3. Global Assembly: Mapping local element DOFs to the global system matrix K
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
                            let val = k_el[(local_node_i * 3 + dof_i, local_node_j * 3 + dof_j)];

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

pub fn compute_element_stiffness(
    project: &Project,
    d_mat: &DMatrix<f64>,
    element: &Element,
) -> Result<DMatrix<f64>, String> {
    let node_ids = element.get_node_ids();
    let num_nodes = node_ids.len();
    let num_dofs = num_nodes * 3;

    // 1. Node Coordinates as nalgebra DMatrix (3 rows x num_nodes columns)
    let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
    for (i, &node_id) in node_ids.iter().enumerate() {
        let coords = project
            .mesh
            .nodes
            .get(&node_id)
            .ok_or(format!("Node {node_id} not found"))?;
        node_coords[(0, i)] = coords.x;
        node_coords[(1, i)] = coords.y;
        node_coords[(2, i)] = coords.z;
    }

    let mut k_el = DMatrix::<f64>::zeros(num_dofs, num_dofs);

    // 2. Numerical Integration Loop over Gauss Points
    for gp in element.integration_points() {
        let (_, dn_local) = element.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);

        // 3. Jacobian Calculation: J = dN_local * node_coords^T
        // Result is a 3x3 matrix (physical space dims x reference space dims)
        let jacobian_dm = &dn_local * node_coords.transpose();

        // Convert to fixed-size 3x3 for faster inversion
        let jacobian = SMatrix::<f64, 3, 3>::from_iterator(jacobian_dm.iter().copied());

        let det_j = jacobian.determinant();
        if det_j.abs() < 1e-14 {
            return Err(format!("Singular element found with nodes: {node_ids:?}"));
        }
        let inv_j = jacobian.try_inverse().unwrap();

        // 4. Transform derivatives to physical space: dN_global = J^-1 * dN_local
        let dn_global = inv_j * dn_local;

        // 5. Build B-Matrix (needs to be updated to return DMatrix)
        let b_mat = build_b_matrix(&dn_global, num_nodes);

        // 6. Core of Weak Form: B^T * D * B
        // We use references (&) to avoid moving/cloning matrices in the loop
        let btdb = b_mat.transpose() * (d_mat * &b_mat);

        // 7. Accumulate: k_el += B^T * D * B * det(J) * w_i
        k_el += btdb * (det_j * gp.weight);
    }

    Ok(k_el)
}

use nalgebra::{Dim, Matrix, storage::Storage};
/// Builds the Strain-Displacement Matrix (B) for a 3D element.
/// Maps nodal displacements to the 6 components of the strain tensor.
pub fn build_b_matrix<R, C, S>(dn_global: &Matrix<f64, R, C, S>, num_nodes: usize) -> DMatrix<f64>
where
    R: Dim,
    C: Dim,
    S: Storage<f64, R, C>,
{
    let num_dofs = num_nodes * 3;
    let mut b = DMatrix::<f64>::zeros(6, num_dofs);

    for i in 0..num_nodes {
        let col_idx = i * 3;

        // nalgebra indexing (row, col) works for any storage type
        let d_dx = dn_global[(0, i)];
        let d_dy = dn_global[(1, i)];
        let d_dz = dn_global[(2, i)];

        // Standard 3D B-Matrix mapping (Voigt notation)
        // Rows: epsilon_xx, epsilon_yy, epsilon_zz, gamma_xy, gamma_yz, gamma_zx
        b[(0, col_idx)] = d_dx;
        b[(1, col_idx + 1)] = d_dy;
        b[(2, col_idx + 2)] = d_dz;

        b[(3, col_idx)] = d_dy;
        b[(3, col_idx + 1)] = d_dx;

        b[(4, col_idx + 1)] = d_dz;
        b[(4, col_idx + 2)] = d_dy;

        b[(5, col_idx)] = d_dz;
        b[(5, col_idx + 2)] = d_dx;
    }
    b
}
