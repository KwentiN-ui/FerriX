#[allow(unused_imports)]
use crate::solver::ids::{NodeId, LoadId, BoundaryConditionId};
use std::{error::Error, time::Duration};

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array2, ArrayView1};
use sprs::{CsMat, TriMat};

use crate::solver::{
    mesh_lib::elements::element::Element,
    project::Project,
    results::{FieldType, NodalResult, StepResult},
    step::boundary_conds::{BoundaryCondition, Load},
};

#[derive(Debug, Clone)]
pub struct StaticStep {
    project: Box<Project>,
}

impl StaticStep {
    pub fn new(project: Box<Project>) -> Self {
        Self { project }
    }

    pub fn compute(&mut self, step_id: usize, loads: &[Load], bcs: &[BoundaryCondition]) -> Result<StepResult, Box<dyn Error>> {
        // 1. Setup
        println!("Constructing global stiffness matrix");

        let num_nodes = self.project.mesh.nodes.len();
        if num_nodes == 0 {
            return Err("Mesh empty or mappings not initialized".into());
        }
        let num_dofs = num_nodes * 3;

        let mut triplet = TriMat::new((num_dofs, num_dofs));

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

        // 3. Assemble stiffness matrix (Element Loop)
        let mut max_diag_val: f64 = 0.0;

        for element in self.project.mesh.elements.values() {
            // --- Get material for this element ---
            let material_index = self
                .project
                .element_materials
                .get(&element.get_id())
                .ok_or(format!(
                    "Element {} has no material assigned.",
                    element.get_id()
                ))?;
            let material = &self.project.materials[*material_index];
            let d_matrix = material.build_elastic_d_matrix();
            // --- End get material ---

            let k_el = self.compute_element_stiffness(&d_matrix, element)?;
            let node_ids = element.get_node_ids();

            for (local_node_i, &global_id_i) in node_ids.iter().enumerate() {
                let global_index_i = self
                    .project
                    .mesh
                    .get_index_for_node_id(global_id_i)
                    .ok_or(format!("Node {global_id_i} not found"))?;

                for (local_node_j, &global_id_j) in node_ids.iter().enumerate() {
                    let global_index_j = self
                        .project
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
            if let Some(idx) = self.project.mesh.get_index_for_node_id(bc.node_id) {
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

        println!(
            "System assembled. K: {}x{}, NNZ: {}. Solving...",
            k_global.rows(),
            k_global.cols(),
            k_global.nnz()
        );

        // Call solver
        let u = solve_cg(&k_global, &f_global, 1e-8, 10000)?;

        // Log result (e.g. displacement norm)
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
                .project
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
fn solve_cg(k: &CsMat<f64>, b: &[f64], tol: f64, max_iter: usize) -> Result<Vec<f64>, String> {
    println!("Conjugate gradient solver is running...");

    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(100));
    spinner.set_style(ProgressStyle::default_spinner());

    let b_len = b.len();
    let mut x = vec![0.0; b_len]; // Initial guess x0 = 0

    // r = b - K * x0  (since x0=0 -> r=b)
    let mut r = b.to_vec();
    let mut p = r.clone();

    // r_old = r^T * r
    let mut rs_old: f64 = r.iter().map(|v| v * v).sum();

    for _ in 0..max_iter {
        spinner.set_message(format!("Residual: {rs_old:.6}"));
        if rs_old.sqrt() < tol {
            spinner.finish();
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
    }
    spinner.finish();
    Err(format!(
        "CG solver did not converge after {max_iter} iterations"
    ))
}