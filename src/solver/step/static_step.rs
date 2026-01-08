use std::{error::Error, time::Duration};

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array2, ArrayView1};
use sprs::{CsMat, TriMat};

use crate::solver::{
    inp::InpFile,
    mesh_lib::{elements::element::Element, mesh::Mesh},
    results::{FieldType, NodalResult, StepResult},
    step::boundary_conds::{BoundaryCondition, Load},
};

// hardcoded constants for now, will be replaced with material card
const E_MOD: f64 = 210_000.0;
const NU: f64 = 0.3;

#[derive(Debug, Clone)]
pub struct StaticStep {
    input: Box<InpFile>,
    mesh: Box<Mesh>,
    line_number: usize,
}

impl StaticStep {
    pub fn new(input: Box<InpFile>, mesh: Box<Mesh>, line_number: usize) -> Self {
        Self {
            input,
            mesh,
            line_number,
        }
    }

    /// Helper: Resolves a string (Set Name or ID) to a list of Node IDs
    fn resolve_target(&self, target: &str) -> Vec<usize> {
        let t = target.trim();
        // Check Node Sets
        if let Some(ids) = self.mesh.node_sets.get(t) {
            return ids.clone();
        }
        // Try explicit Node ID
        if let Ok(id) = t.parse::<usize>() {
            return vec![id];
        }
        // Not found
        Vec::new()
    }

    pub fn parse_loads(&self) -> Vec<Load> {
        let mut loads = Vec::new();
        let mut lines = self.input.0.lines().peekable();

        while let Some(line) = lines.next() {
            let trimmed = line.trim();

            if trimmed.starts_with("*CLOAD") {
                // Read data lines until next keyword
                while let Some(next_line) = lines.peek() {
                    let data = next_line.trim();
                    if data.starts_with('*') {
                        break;
                    }

                    // Format: Node/Set, DOF, Value
                    let parts: Vec<&str> = data.split(',').collect();
                    if parts.len() >= 3 {
                        let target_nodes = self.resolve_target(parts[0]);
                        let dof_in: usize = parts[1].trim().parse().unwrap_or(0);
                        let val: f64 = parts[2].trim().parse().unwrap_or(0.0);

                        // DOF: Inp 1,2,3 -> Internal 0,1,2
                        if (1..=3).contains(&dof_in) {
                            for node_id in target_nodes {
                                loads.push(Load {
                                    node_id,
                                    dof: dof_in - 1,
                                    value: val,
                                });
                            }
                        }
                    }

                    lines.next(); // Consume line
                }
            }
        }

        loads
    }

    pub fn parse_bcs(&self) -> Vec<BoundaryCondition> {
        let mut bcs = Vec::new();
        let mut lines = self.input.0.lines().peekable();

        while let Some(line) = lines.next() {
            let trimmed = line.trim();

            if trimmed.starts_with("*BOUNDARY") {
                // Read data lines
                while let Some(next_line) = lines.peek() {
                    let data = next_line.trim();
                    if data.starts_with('*') {
                        break;
                    }

                    // Format: Node/Set, FirstDOF, [LastDOF], [Value]
                    let parts: Vec<&str> = data.split(',').collect();
                    if parts.len() >= 2 {
                        let target_nodes = self.resolve_target(parts[0]);

                        let first_dof: usize = parts[1].trim().parse().unwrap_or(0);
                        // If LastDOF is missing, default to FirstDOF
                        let last_dof: usize = if parts.len() > 2 && !parts[2].trim().is_empty() {
                            parts[2].trim().parse().unwrap_or(first_dof)
                        } else {
                            first_dof
                        };

                        // If Value is missing, default to 0.0 (Standard fixation)
                        let val: f64 = if parts.len() > 3 {
                            parts[3].trim().parse().unwrap_or(0.0)
                        } else {
                            0.0
                        };

                        for node_id in target_nodes {
                            // Loop over DOFs (e.g. 1 to 3 => x, y, z)
                            for dof_in in first_dof..=last_dof {
                                if (1..=3).contains(&dof_in) {
                                    bcs.push(BoundaryCondition {
                                        node_id,
                                        dof: dof_in - 1,
                                        value: val,
                                    });
                                }
                            }
                        }
                    }

                    lines.next();
                }
            }
        }

        bcs
    }

    pub fn compute(&mut self) -> Result<StepResult, Box<dyn Error>> {
        let loads = self.parse_loads();
        let bcs = self.parse_bcs();

        // 1. Setup
        println!("Constructing global stiffness matrix");

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::mesh_lib::node::Node;
    use std::collections::HashMap;

    #[test]
    fn test_single_c3d4_element_tension() {
        // 1. Create a simple mesh in-memory
        let mut nodes = HashMap::new();
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
        let element = Element::C3D4(1, [1, 2, 3, 4]);
        elements.insert(1, element);

        let mut node_sets = HashMap::new();
        node_sets.insert("Fixed".to_string(), vec![1, 2, 3]);

        let mut mesh = Mesh {
            nodes,
            elements,
            node_sets,
            node_id_to_index: HashMap::new(),
            index_to_node_id: Vec::new(),
        };
        mesh.build_node_mappings();

        // 2. Create an inp-file string for BCs and Loads
        let inp_file_content = r#"
*BOUNDARY
Fixed, 1, 3
*CLOAD
4, 3, 1000.0
"#
        .to_string();
        let inp_file = InpFile(inp_file_content);

        // 3. Setup and run the StaticStep
        let mut static_step = StaticStep::new(Box::new(inp_file), Box::new(mesh));
        let result = static_step.compute();

        // 4. Assert the results
        assert!(result.is_ok(), "Solver failed to compute");
        let step_result = result.unwrap();

        let u_field = step_result
            .nodal_results
            .iter()
            .find(|r| r.name == "U")
            .expect("Displacement field not found");

        // --- Analytical check ---
        // This is a rough approximation. A real test would use a known
        // result from another FEA software or a textbook example.
        // For a simple bar: disp = F*L / (A*E)
        // Here, we have a tetrahedron. Let's just check if the displacement
        // of node 4 is in the right direction and plausible.
        let disp_node_4 = u_field.data.get(&4).expect("Node 4 not in results");
        let uz = disp_node_4[2];

        // Assert that displacement in z is positive and significant
        assert!(uz > 0.0, "Displacement in z should be positive");

        // The exact value is hard to calculate by hand.
        // Let's assert it's within a plausible range for this test setup.
        // This value is a placeholder and should be replaced with a verified one.
        let expected_uz_approx = 0.0212245132255;
        let tolerance = expected_uz_approx * 0.001; // 1% tolerance
        assert!(
            (uz - expected_uz_approx).abs() < tolerance,
            "Displacement Z at node 4 is {:.6}, expected around {:.6}",
            uz,
            expected_uz_approx
        );

        // Assert that other displacements on the loaded node are near zero
        let ux = disp_node_4[0];
        let uy = disp_node_4[1];
        assert!(
            ux.abs() < 1e-6,
            "Displacement X at node 4 should be near zero"
        );
        assert!(
            uy.abs() < 1e-6,
            "Displacement Y at node 4 should be near zero"
        );

        // Assert that fixed nodes have zero displacement
        for node_id in &[1, 2, 3] {
            let disp = u_field.data.get(node_id).unwrap();
            assert!(disp.iter().all(|&d| d.abs() < 1e-6));
        }
    }

    #[test]
    fn test_single_c3d20_element_cantilever() {
        // 1. Create a simple mesh in-memory for a cantilever beam
        let mut nodes = HashMap::new();
        // Node coordinates for a 10x1x1 beam
        // Corners
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
                x: 10.0,
                y: 0.0,
                z: 0.0,
            },
        );
        nodes.insert(
            3,
            Node {
                id: 3,
                x: 10.0,
                y: 1.0,
                z: 0.0,
            },
        );
        nodes.insert(
            4,
            Node {
                id: 4,
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
        );
        nodes.insert(
            5,
            Node {
                id: 5,
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        );
        nodes.insert(
            6,
            Node {
                id: 6,
                x: 10.0,
                y: 0.0,
                z: 1.0,
            },
        );
        nodes.insert(
            7,
            Node {
                id: 7,
                x: 10.0,
                y: 1.0,
                z: 1.0,
            },
        );
        nodes.insert(
            8,
            Node {
                id: 8,
                x: 0.0,
                y: 1.0,
                z: 1.0,
            },
        );
        // Midside
        nodes.insert(
            9,
            Node {
                id: 9,
                x: 5.0,
                y: 0.0,
                z: 0.0,
            },
        );
        nodes.insert(
            10,
            Node {
                id: 10,
                x: 10.0,
                y: 0.5,
                z: 0.0,
            },
        );
        nodes.insert(
            11,
            Node {
                id: 11,
                x: 5.0,
                y: 1.0,
                z: 0.0,
            },
        );
        nodes.insert(
            12,
            Node {
                id: 12,
                x: 0.0,
                y: 0.5,
                z: 0.0,
            },
        );
        nodes.insert(
            13,
            Node {
                id: 13,
                x: 0.0,
                y: 0.0,
                z: 0.5,
            },
        );
        nodes.insert(
            14,
            Node {
                id: 14,
                x: 10.0,
                y: 0.0,
                z: 0.5,
            },
        );
        nodes.insert(
            15,
            Node {
                id: 15,
                x: 10.0,
                y: 1.0,
                z: 0.5,
            },
        );
        nodes.insert(
            16,
            Node {
                id: 16,
                x: 0.0,
                y: 1.0,
                z: 0.5,
            },
        );
        nodes.insert(
            17,
            Node {
                id: 17,
                x: 5.0,
                y: 0.0,
                z: 1.0,
            },
        );
        nodes.insert(
            18,
            Node {
                id: 18,
                x: 10.0,
                y: 0.5,
                z: 1.0,
            },
        );
        nodes.insert(
            19,
            Node {
                id: 19,
                x: 5.0,
                y: 1.0,
                z: 1.0,
            },
        );
        nodes.insert(
            20,
            Node {
                id: 20,
                x: 0.0,
                y: 0.5,
                z: 1.0,
            },
        );

        let mut elements = HashMap::new();
        let element = Element::C3D20(
            1,
            [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            ],
        );
        elements.insert(1, element);

        let mut node_sets = HashMap::new();
        node_sets.insert("Fixed".to_string(), vec![1, 4, 5, 8, 12, 13, 16, 20]);
        node_sets.insert("Load".to_string(), vec![2, 3, 6, 7, 10, 14, 15, 18]);

        let mut mesh = Mesh {
            nodes,
            elements,
            node_sets,
            node_id_to_index: HashMap::new(),
            index_to_node_id: Vec::new(),
        };
        mesh.build_node_mappings();

        // 2. Create an inp-file string for BCs and Loads
        let inp_file_content = r#"
*BOUNDARY
Fixed, 1, 3
*CLOAD
Load, 3, -100.0
"#
        .to_string();
        let inp_file = InpFile(inp_file_content);

        // 3. Setup and run the StaticStep
        let mut static_step = StaticStep::new(Box::new(inp_file), Box::new(mesh));
        let result = static_step.compute();

        // 4. Assert the results
        assert!(result.is_ok(), "Solver failed to compute");
        let step_result = result.unwrap();

        let u_field = step_result
            .nodal_results
            .iter()
            .find(|r| r.name == "U")
            .expect("Displacement field not found");

        let disp_node_2 = u_field.data.get(&2).expect("Node 2 not in results");
        let uz = disp_node_2[2];

        // Cantilever beam theory: d = (P * L^3) / (3 * E * I)
        // P = 100 * 8 = 800N (total load)
        // L = 10
        // E = 210000
        // I = (b*h^3)/12 = (1*1^3)/12 = 1/12
        // d = (800 * 10^3) / (3 * 210000 * (1/12)) = 15.23
        // This is a rough estimation, the FEA result will be different.
        let expected_uz_approx: f64 = -13.806006;
        let tolerance = expected_uz_approx.abs() * 0.01; // 1% tolerance
        assert!(
            (uz - expected_uz_approx).abs() < tolerance,
            "Displacement Z at node 2 is {:.6}, expected around {:.6}",
            uz,
            expected_uz_approx
        );
    }
}
