use std::{
    collections::HashMap,
    error::Error,
    sync::{Arc, mpsc::Sender},
};

use ndarray::{Array2, ArrayView1};
use sprs::{CsMat, TriMat};

use crate::{
    solver::{
        inp::InpFile,
        mesh_lib::{elements::element::Element, mesh::Mesh},
        results::{FieldType, NodalResult, StepResult},
        step::boundary_conds::{BoundaryCondition, Load},
    },
    tui::app::AppEvent,
};

// hardcoded constants for now, will be replaced with material card
const E_MOD: f64 = 210_000.0;
const NU: f64 = 0.3;

#[derive(Debug, Clone)]
pub struct StaticStep {
    input: Arc<InpFile>,
    mesh: Arc<Mesh>,
}

impl StaticStep {
    pub fn new(input: Arc<InpFile>, mesh: Arc<Mesh>) -> Self {
        Self { input, mesh }
    }

    pub fn parse_loads(&self) -> Vec<Load> {
        // TODO: this is a test!
        let min_id = self.mesh.nodes.keys().min().copied().unwrap_or(1);
        let max_id = self.mesh.nodes.keys().max().copied().unwrap_or(1);

        vec![Load {
            node_id: max_id,
            dof: 0,
            value: 1000.0,
        }]
    }
    pub fn parse_bcs(&self) -> Vec<BoundaryCondition> {
        // TODO: this is a test!
        let fixed_nodes: Vec<usize> = self
            .mesh
            .nodes
            .values()
            .filter(|n| n.x.abs() < 1e-5)
            .map(|n| n.id)
            .collect();

        let mut bcs = Vec::new();
        for id in fixed_nodes {
            // "Festgeschweißt": Alles 0
            bcs.push(BoundaryCondition {
                node_id: id,
                dof: 0,
                value: 0.0,
            });
            bcs.push(BoundaryCondition {
                node_id: id,
                dof: 1,
                value: 0.0,
            });
            bcs.push(BoundaryCondition {
                node_id: id,
                dof: 2,
                value: 0.0,
            });
        }
        bcs
    }

    pub fn compute(&mut self, tx: &Sender<AppEvent>) -> Result<StepResult, Box<dyn Error>> {
        let loads = self.parse_loads();
        let bcs = self.parse_bcs();
        // 1. Setup
        let _ = tx.send(AppEvent::SolverLog(
            "Constructing global stiffness matrix".to_string(),
        ));
        let num_nodes = self.mesh.index_to_node_id.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }
        let num_dofs = num_nodes * 3;

        let mut triplet = TriMat::new((num_dofs, num_dofs));
        let d_matrix = build_elastic_d_matrix(E_MOD, NU);

        // Init Force Vector F
        let mut f_global = vec![0.0; num_dofs];

        // 2. Add loads to F
        for load in loads {
            if let Some(idx) = self.mesh.get_index_for_node_id(load.node_id) {
                let global_dof = idx * 3 + load.dof;
                if global_dof < num_dofs {
                    f_global[global_dof] += load.value;
                }
            } else {
                eprintln!("Warning: Load on unknown node {}", load.node_id);
            }
        }

        // 3. Assemble stiffness matrix (Element Loop)
        let mut max_diag_val: f64 = 0.0;

        for element in self.mesh.elements.values() {
            let k_el = self.compute_element_stiffness(&d_matrix, element)?;
            let node_ids = element.get_node_ids();

            for (local_node_i, &global_id_i) in node_ids.iter().enumerate() {
                let global_index_i = self
                    .mesh
                    .get_index_for_node_id(global_id_i)
                    .ok_or(format!("Node {global_id_i} not found"))?;

                for (local_node_j, &global_id_j) in node_ids.iter().enumerate() {
                    let global_index_j = self
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

                                // Track max diagonal for penalty factor estimation
                                if row == col {
                                    max_diag_val = max_diag_val.max(val.abs());
                                }
                            }
                        }
                    }
                }
            }
        }

        // 4. Apply boundary conditions (Penalty Method)
        // Add a large value to the diagonal in the triplet matrix.
        // This is very efficient for sparse matrices.
        if max_diag_val == 0.0 {
            max_diag_val = 1.0;
        } // Fallback
        let penalty = max_diag_val * 1.0e6;

        for bc in bcs {
            if let Some(idx) = self.mesh.get_index_for_node_id(bc.node_id) {
                let global_dof = idx * 3 + bc.dof;

                if global_dof < num_dofs {
                    // K_ii += alpha
                    triplet.add_triplet(global_dof, global_dof, penalty);

                    // F_i += alpha * u_bc
                    f_global[global_dof] += penalty * bc.value;
                }
            }
        }

        // 5. Conversion & Solving
        let k_global: CsMat<f64> = triplet.to_csr();

        let _ = tx.send(AppEvent::SolverLog(format!(
            "System assembled. K: {}x{}, NNZ: {}. Solving...",
            k_global.rows(),
            k_global.cols(),
            k_global.nnz()
        )));

        // Call solver
        let u = solve_cg(&k_global, &f_global, 1e-8, 10000, tx)?;

        // Log result (e.g. displacement norm)
        let u_norm: f64 = u.iter().map(|x| x * x).sum::<f64>().sqrt();
        let _ = tx.send(AppEvent::SolverLog(format!(
            "Solution converged. Displacement Norm: {u_norm:.4e}"
        )));

        let mut displacement_field =
            NodalResult::new("Static_Displacement", FieldType::Displacement);

        for (matrix_idx, &node_id) in self.mesh.index_to_node_id.iter().enumerate() {
            let idx = matrix_idx * 3;
            if idx + 2 < u.len() {
                let dx = u[idx];
                let dy = u[idx + 1];
                let dz = u[idx + 2];
                displacement_field.insert(node_id, vec![dx, dy, dz]);
            }
        }

        let mut step_res = StepResult::new(1, "Static Step", 1.0);
        step_res.nodal_results.push(displacement_field);

        Ok(step_res)
    }

    fn compute_element_stiffness(
        &self,
        d_mat: &Array2<f64>,
        element: &Element,
    ) -> Result<Array2<f64>, String> {
        let node_ids = element.get_node_ids();
        let num_nodes = node_ids.len();
        let num_dofs = num_nodes * 3;

        // get node coords
        let mut node_coords = Array2::<f64>::zeros((3, num_nodes));
        for (i, &node_id) in node_ids.iter().enumerate() {
            let coords = self
                .mesh
                .nodes
                .get(&node_id)
                .ok_or(format!("Node {node_id} not found"))?;
            node_coords[[0, i]] = coords.x;
            node_coords[[1, i]] = coords.y;
            node_coords[[2, i]] = coords.z;
        }

        let mut k_el = Array2::<f64>::zeros((num_dofs, num_dofs));

        // integration loop
        for gp in element.integration_points() {
            // shape function & local derivative (3 x N)
            let (_, dn_local) = element.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);

            // Jacobi-Matrix: J = dN_local * NodeCoords^T
            let jacobian = dn_local.dot(&node_coords.t());

            // invert & determinant
            let (det_j, inv_j) = invert_jacobian_3x3(&jacobian)
                .map_err(|()| format!("Singular element found with nodes: {node_ids:?}"))?;

            // global diff: dN_global = J^-1 * dN_local
            let dn_global = inv_j.dot(&dn_local);

            // B-Matrix
            let b_mat = build_b_matrix(&dn_global, num_nodes);

            // sum stiffness: K += B^T * D * B * detJ * weight
            let db = d_mat.dot(&b_mat);
            let btdb = b_mat.t().dot(&db);

            // k_el = k_el + scaled_matrix
            k_el.scaled_add(det_j * gp.weight, &btdb);
        }

        Ok(k_el)
    }
}

/// Constructs the B-Matrix (6 x 3*Nodes)
/// using Voigt notation: xx, yy, zz, xy, yz, zx
fn build_b_matrix(dn_global: &Array2<f64>, num_nodes: usize) -> Array2<f64> {
    let num_dofs = num_nodes * 3;
    let mut b = Array2::<f64>::zeros((6, num_dofs));

    for i in 0..num_nodes {
        let col_idx = i * 3;
        let d_dx = dn_global[[0, i]];
        let d_dy = dn_global[[1, i]];
        let d_dz = dn_global[[2, i]];

        // Row 0: epsilon_xx -> dN/dx at u_x
        b[[0, col_idx]] = d_dx;

        // Row 1: epsilon_yy -> dN/dy at u_y
        b[[1, col_idx + 1]] = d_dy;

        // Row 2: epsilon_zz -> dN/dz at u_z
        b[[2, col_idx + 2]] = d_dz;

        // Row 3: gamma_xy -> dN/dy at u_x + dN/dx at u_y
        b[[3, col_idx]] = d_dy;
        b[[3, col_idx + 1]] = d_dx;

        // Row 4: gamma_yz -> dN/dz at u_y + dN/dy at u_z
        b[[4, col_idx + 1]] = d_dz;
        b[[4, col_idx + 2]] = d_dy;

        // Row 5: gamma_zx -> dN/dz at u_x + dN/dx at u_z
        b[[5, col_idx]] = d_dz;
        b[[5, col_idx + 2]] = d_dx;
    }
    b
}

/// constructs linear elastic materialmatrix D (6x6 for 3D)
/// Voigt-Notation: xx, yy, zz, xy, yz, zx
fn build_elastic_d_matrix(e: f64, nu: f64) -> Array2<f64> {
    let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let c1 = 1.0 - nu;
    let c2 = nu;
    let c3 = (1.0 - 2.0 * nu) / 2.0;

    let data = vec![
        c1, c2, c2, 0., 0., 0., c2, c1, c2, 0., 0., 0., c2, c2, c1, 0., 0., 0., 0., 0., 0., c3, 0.,
        0., 0., 0., 0., 0., c3, 0., 0., 0., 0., 0., 0., c3,
    ];

    Array2::from_shape_vec((6, 6), data).expect("Matrix shape error") * factor
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

/// Simple Conjugate Gradient Solver for sparse symmetric positive definite systems.
/// Solves K * x = b
#[allow(clippy::many_single_char_names)]
fn solve_cg(
    k: &CsMat<f64>,
    b: &[f64],
    tol: f64,
    max_iter: usize,
    tx: &Sender<AppEvent>,
) -> Result<Vec<f64>, String> {
    let b_len = b.len();
    let mut x = vec![0.0; b_len]; // Initial guess x0 = 0

    // r = b - K * x0  (since x0=0 -> r=b)
    let mut r = b.to_vec();
    let mut p = r.clone();

    // r_old = r^T * r
    let mut rs_old: f64 = r.iter().map(|v| v * v).sum();

    for i in 0..max_iter {
        let _ = tx.send(AppEvent::SolverLog(format!("Iteration {i}")));
        if rs_old.sqrt() < tol {
            return Ok(x);
        }

        // Ap = K * p
        let p_view = ArrayView1::from(&p);
        let ap = (k * &p_view).to_vec();

        // alpha = rs_old / (p^T * Ap)
        let p_dot_ap: f64 = p.iter().zip(&ap).map(|(pi, api)| pi * api).sum();
        if p_dot_ap.abs() < 1e-15 {
            return Err("CG breakdown: denominator zero (matrix singular?)".to_string());
        }
        let alpha = rs_old / p_dot_ap;

        // x = x + alpha * p
        // r = r - alpha * Ap
        for j in 0..b_len {
            x[j] += alpha * p[j];
            r[j] -= alpha * ap[j];
        }

        let rs_new: f64 = r.iter().map(|v| v * v).sum();

        // beta = rs_new / rs_old
        let beta = rs_new / rs_old;

        // p = r + beta * p
        for j in 0..b_len {
            p[j] = r[j] + beta * p[j];
        }

        rs_old = rs_new;
        let _ = tx.send(AppEvent::SolverLog(format!("Residual: {rs_new}")));
    }

    Err(format!(
        "CG solver did not converge after {max_iter} iterations"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::mesh_lib::elements::element::{Element, ElementType};
    use crate::solver::mesh_lib::node::Node;
    use std::collections::HashMap;

    // Helper to create a dummy step with a single element
    fn create_single_element_step() -> (StaticStep, usize) {
        let mut nodes = HashMap::new();
        // Unit Tetrahedral: Origin + Unit vectors
        nodes.insert(
            1,
            Node {
                id: 1,
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        );
        nodes.insert(
            2,
            Node {
                id: 2,
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
        );
        nodes.insert(
            3,
            Node {
                id: 3,
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
        );
        nodes.insert(
            4,
            Node {
                id: 4,
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        );

        let mut elements = HashMap::new();
        let el_id = 1;
        // C3D4 connecting nodes 1, 2, 3, 4
        let element = Element::C3D4(el_id, [1, 2, 3, 4]);
        elements.insert(el_id, element);

        let mut mesh = Mesh::new();
        mesh.nodes = nodes;
        mesh.elements = elements;
        mesh.build_node_mappings(); // Crucial: Initialize mappings

        // Dummy input file
        let inp_file = InpFile(String::new());

        let step = StaticStep::new(Arc::new(inp_file), Arc::new(mesh));
        (step, el_id)
    }

    #[test]
    fn test_c3d4_stiffness_properties() {
        let (step, el_id) = create_single_element_step();

        // Standard steel parameters
        let e_mod = 210000.0;
        let nu = 0.3;
        let d_matrix = build_elastic_d_matrix(e_mod, nu);

        let element = step.mesh.elements.get(&el_id).unwrap();

        // Calculate Stiffness Matrix
        let k_el = step
            .compute_element_stiffness(&d_matrix, element)
            .expect("Computation failed");

        // 1. Check Dimensions (4 Nodes * 3 DOFs = 12x12)
        assert_eq!(
            k_el.shape(),
            &[12, 12],
            "Stiffness matrix has wrong dimensions"
        );

        // 2. Check Symmetry (K_ij == K_ji)
        for i in 0..12 {
            for j in i + 1..12 {
                let val_ij = k_el[[i, j]];
                let val_ji = k_el[[j, i]];
                assert!(
                    (val_ij - val_ji).abs() < 1e-9,
                    "Matrix not symmetric at ({}, {}): {} != {}",
                    i,
                    j,
                    val_ij,
                    val_ji
                );
            }
        }

        // 3. Check Rigid Body Translation
        // Sum of rows must be zero (equilibrium of forces for rigid move)
        for i in 0..12 {
            let row_sum: f64 = k_el.row(i).sum();
            assert!(
                row_sum.abs() < 1e-9,
                "Row {} sum is not zero ({}), rigid body motion fails",
                i,
                row_sum
            );
        }

        // 4. Check Positive Diagonal
        // Diagonal elements represent stiffness against direct displacement, usually > 0
        for i in 0..12 {
            assert!(k_el[[i, i]] > 0.0, "Diagonal element {} is not positive", i);
        }
    }
}
