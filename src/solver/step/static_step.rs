use std::{
    error::Error,
    sync::{Arc, mpsc::Sender},
};

use ndarray::Array2;
use sprs::{CsMat, TriMat};

use crate::{
    solver::{
        inp::InpFile,
        mesh_lib::{elements::element::Element, mesh::Mesh},
    },
    tui::app::AppEvent,
};

// Konstanten für Standardstahl (erstmal hardcoded)
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

    pub fn compute(&mut self, tx: &Sender<AppEvent>) -> Result<(), Box<dyn Error>> {
        let _ = tx.send(AppEvent::SolverLog(
            "Assembling global stiffness matrix".to_string(),
        ));
        // 1. DOF Management
        let num_nodes = self.mesh.nodes.len();
        let num_dofs = num_nodes * 3; // 3 DOFs pro Knoten (u_x, u_y, u_z)

        // Triplet Matrix initialisieren (Row, Col, Value)
        let mut triplet = TriMat::new((num_dofs, num_dofs));

        // 2. Materialmatrix D berechnen (ndarray)
        let d_matrix = build_elastic_d_matrix(E_MOD, NU);

        // 3. Element-Loop
        // Fix: Da .elements eine Map ist, bekommen wir (id, element)
        for element in self.mesh.elements.values() {
            // Lokale Steifigkeitsmatrix K_el berechnen
            // Fix: Das '?' am Ende behandelt das Result, falls das Element singulär ist
            let k_el = self.compute_element_stiffness(&d_matrix, element)?;

            // Node IDs für das Mapping holen
            // Fix: Aufruf der Methode get_node_ids() statt Feldzugriff
            let node_ids = element.get_node_ids();

            // 4. Assemblierung: Local -> Global
            for (local_node_i, &global_node_i) in node_ids.iter().enumerate() {
                for (local_node_j, &global_node_j) in node_ids.iter().enumerate() {
                    // Für jeden Knoten 3 DOFs durchgehen
                    for dof_i in 0..3 {
                        for dof_j in 0..3 {
                            // Mapping auf globale Gleichungsnummern
                            let global_row = global_node_i * 3 + dof_i;
                            let global_col = global_node_j * 3 + dof_j;

                            // Wert aus lokalem k_el holen
                            let val = k_el[[local_node_i * 3 + dof_i, local_node_j * 3 + dof_j]];

                            if val.abs() > 1e-12 {
                                triplet.add_triplet(global_row, global_col, val);
                            }
                        }
                    }
                }
            }
        }

        // 5. Konvertierung in CSR (Compressed Sparse Row) für den Solver
        let k_global: CsMat<f64> = triplet.to_csr();

        println!(
            "Assemblierung fertig. K: ({}, {}) mit {} Einträgen",
            k_global.rows(),
            k_global.cols(),
            k_global.nnz()
        );

        Ok(())
    }

    fn compute_element_stiffness(
        &self,
        d_mat: &Array2<f64>,
        element: &Element,
    ) -> Result<Array2<f64>, String> {
        let node_ids = element.get_node_ids();
        let num_nodes = node_ids.len();
        let num_dofs = num_nodes * 3;

        // 1. Knoten-Koordinaten holen
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

        // 2. Integration Loop
        for gp in element.integration_points() {
            // Formfunktionen & Lokale Ableitungen (3 x N)
            let (_, dn_local) = element.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);

            // Jacobi-Matrix: J = dN_local * NodeCoords^T
            // Achtung: node_coords transponieren für korrekte Multiplikation
            let jacobian = dn_local.dot(&node_coords.t());

            // Invertierung & Determinante
            let (det_j, inv_j) = invert_jacobian_3x3(&jacobian)
                .map_err(|()| format!("Singular element found with nodes: {node_ids:?}"))?;

            // Globale Ableitungen: dN_global = J^-1 * dN_local
            let dn_global = inv_j.dot(&dn_local);

            // B-Matrix aufstellen
            let b_mat = build_b_matrix(&dn_global, num_nodes);

            // Steifigkeit aufaddieren: K += B^T * D * B * detJ * weight
            let db = d_mat.dot(&b_mat);
            let btdb = b_mat.t().dot(&db);

            // k_el = k_el + scaled_matrix
            k_el.scaled_add(det_j * gp.weight, &btdb);
        }

        Ok(k_el)
    }
}

/// Erstellt die B-Matrix (6 x 3*Nodes)
/// Reihenfolge der Strains (Voigt): xx, yy, zz, xy, yz, zx
fn build_b_matrix(dn_global: &Array2<f64>, num_nodes: usize) -> Array2<f64> {
    let num_dofs = num_nodes * 3;
    let mut b = Array2::<f64>::zeros((6, num_dofs));

    for i in 0..num_nodes {
        let col_idx = i * 3;
        let d_dx = dn_global[[0, i]];
        let d_dy = dn_global[[1, i]];
        let d_dz = dn_global[[2, i]];

        // Zeile 0: epsilon_xx -> dN/dx an u_x
        b[[0, col_idx]] = d_dx;

        // Zeile 1: epsilon_yy -> dN/dy an u_y
        b[[1, col_idx + 1]] = d_dy;

        // Zeile 2: epsilon_zz -> dN/dz an u_z
        b[[2, col_idx + 2]] = d_dz;

        // Zeile 3: gamma_xy -> dN/dy an u_x + dN/dx an u_y
        b[[3, col_idx]] = d_dy;
        b[[3, col_idx + 1]] = d_dx;

        // Zeile 4: gamma_yz -> dN/dz an u_y + dN/dy an u_z
        b[[4, col_idx + 1]] = d_dz;
        b[[4, col_idx + 2]] = d_dy;

        // Zeile 5: gamma_zx -> dN/dz an u_x + dN/dx an u_z
        b[[5, col_idx]] = d_dz;
        b[[5, col_idx + 2]] = d_dx;
    }
    b
}

/// Erstellt die linearelastische Materialmatrix D (6x6 für 3D)
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

/// Helper: Invertiert eine 3x3 Matrix und gibt Determinante zurück.
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
