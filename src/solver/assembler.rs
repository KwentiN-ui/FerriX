use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::project::Project;
use nalgebra::{DMatrix, SMatrix};
use sprs::TriMat;
use std::ops::AddAssign;

pub struct Assembler;

impl Assembler {
    pub fn assemble(project: &Project, is_symmetric: bool) -> Result<(TriMat<f64>, f64), String> {
        let num_nodes = project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }

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

            let d_matrix = material.build_elastic_d_matrix();
            let k_el = compute_element_stiffness(project, &d_matrix, element, is_symmetric)?;
            let node_ids = element.get_node_ids();
            let num_nodes_el = node_ids.len();

            for i in 0..num_nodes_el {
                let global_index_i = project
                    .mesh
                    .get_index_for_node_id(node_ids[i])
                    .ok_or(format!("Node {} not found", node_ids[i]))?;

                // If symmetric, only iterate from i to avoid redundant calculations
                let start_j = if is_symmetric { i } else { 0 };

                for j in start_j..num_nodes_el {
                    let global_index_j = project
                        .mesh
                        .get_index_for_node_id(node_ids[j])
                        .ok_or(format!("Node {} not found", node_ids[j]))?;

                    for dof_i in 0..3 {
                        // If symmetric and on the same node, only iterate from dof_i
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
                                    // Add the transposed entry for symmetric matrices
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

pub fn compute_element_stiffness(
    project: &Project,
    d_mat: &DMatrix<f64>,
    element: &Element,
    is_symmetric: bool,
) -> Result<DMatrix<f64>, String> {
    let node_ids = element.get_node_ids();
    let num_nodes = node_ids.len();
    let num_dofs = num_nodes * 3;

    let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
    for (i, &node_id) in node_ids.iter().enumerate() {
        let coords = project.mesh.nodes.get(&node_id).ok_or("Node not found")?;
        node_coords[(0, i)] = coords.x;
        node_coords[(1, i)] = coords.y;
        node_coords[(2, i)] = coords.z;
    }

    let mut k_el = DMatrix::<f64>::zeros(num_dofs, num_dofs);

    for gp in element.integration_points() {
        let (_, dn_local) = element.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
        let jacobian = &dn_local * node_coords.transpose();
        let det_j = jacobian.determinant();
        let inv_j = jacobian.try_inverse().ok_or("Singular Jacobian")?;

        let dn_global = inv_j * dn_local;
        let weight = det_j.abs() * gp.weight;

        for i in 0..num_nodes {
            let bi = build_b_block(&dn_global, i);

            // If symmetric: only j from i..num_nodes
            // If asymmetric: j from 0..num_nodes
            let start_j = if is_symmetric { i } else { 0 };

            for j in start_j..num_nodes {
                let bj = build_b_block(&dn_global, j);
                let k_block = bi.transpose() * (d_mat * bj) * weight;

                k_el.fixed_view_mut::<3, 3>(i * 3, j * 3)
                    .add_assign(&k_block);

                if is_symmetric && i != j {
                    // Exploit $k_{ji} = k_{ij}^T$
                    k_el.fixed_view_mut::<3, 3>(j * 3, i * 3)
                        .add_assign(&k_block.transpose());
                } else if !is_symmetric && i != j {
                    // In asymmetric case, the loop over j starts at 0,
                    // so k_ji will be calculated naturally in another iteration.
                }
            }
        }
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

/// Helper to get the 6x3 sub-block of B for a single node (Voigt notation)
fn build_b_block(dn_global: &DMatrix<f64>, node_idx: usize) -> SMatrix<f64, 6, 3> {
    let mut b = SMatrix::<f64, 6, 3>::zeros();
    let dx = dn_global[(0, node_idx)];
    let dy = dn_global[(1, node_idx)];
    let dz = dn_global[(2, node_idx)];

    // Strains: [eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx]^T
    b[(0, 0)] = dx;
    b[(1, 1)] = dy;
    b[(2, 2)] = dz;

    b[(3, 0)] = dy;
    b[(3, 1)] = dx;
    b[(4, 1)] = dz;
    b[(4, 2)] = dy;
    b[(5, 0)] = dz;
    b[(5, 2)] = dx;
    b
}
