use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::project::Project;
use nalgebra::{DMatrix, DVector, SMatrix};
use std::ops::AddAssign;
use std::str::FromStr;

use crate::solver::mesh_lib::elements::c3d4::{c3d4_gauss, shape_func_c3d4};
use strum_macros::{EnumDiscriminants, EnumString};

#[derive(EnumDiscriminants, Debug, Clone)]
#[strum_discriminants(derive(Hash, EnumString))]
#[strum_discriminants(name(ElementType))]
pub enum Element {
    C3D4(ElementId, [NodeId; 4]),
}

impl Element {
    /// Parses the element type from a string.
    ///
    /// # Errors
    /// Returns an error if the element type is not found.
    pub fn parse_type_str_from_line(line: &str) -> Result<String> {
        Ok(line
            .split(',')
            .map(str::trim)
            .nth(1)
            .ok_or_else(|| FerrixError::ParseError {
                line: 0,
                message: "Invalid element definition".into(),
            })?
            .split('=')
            .next_back()
            .ok_or_else(|| FerrixError::ParseError {
                line: 0,
                message: "Invalid element definition".into(),
            })?
            .to_string())
    }

    /// Parses an element from a string line.
    ///
    /// # Errors
    /// Returns an error if the input line is malformed or an unknown element type is encountered.
    pub fn parse_line(type_name: &str, line: &str) -> Result<Self> {
        let nums: Vec<usize> = line
            .split(',')
            .map(|s| {
                s.trim().parse().map_err(|_| FerrixError::ParseError {
                    line: 0,
                    message: "Integer conversion failed".into(),
                })
            })
            .collect::<Result<Vec<usize>>>()?;

        let (&id, nodes_usize) = nums.split_first().ok_or_else(|| FerrixError::ParseError {
            line: 0,
            message: "Line empty".into(),
        })?;
        let id = ElementId(id);
        let nodes: Vec<NodeId> = nodes_usize.iter().map(|&n| NodeId(n)).collect();

        let elem_type = ElementType::from_str(type_name).map_err(|_| FerrixError::ParseError {
            line: 0,
            message: format!("Unknown element definition: {type_name}"),
        })?;

        match elem_type {
            ElementType::C3D4 => {
                let nodes_arr: [NodeId; 4] =
                    nodes.try_into().map_err(|_| FerrixError::ParseError {
                        line: 0,
                        message: "Wrong node count for C3D4".into(),
                    })?;
                Ok(Element::C3D4(id, nodes_arr))
            }
        }
    }

    #[must_use]
    pub fn get_id(&self) -> ElementId {
        match self {
            Element::C3D4(id, _) => *id,
        }
    }

    #[must_use]
    pub fn get_node_ids(&self) -> &[NodeId] {
        match self {
            Element::C3D4(_, n) => n,
        }
    }

    #[must_use]
    pub fn integration_points(&self) -> Vec<GaussPoint> {
        match self {
            Element::C3D4(..) => c3d4_gauss(),
        }
    }

    #[must_use]
    pub fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        match self {
            Element::C3D4(..) => {
                // Get static shape functions and derivatives
                let (n, dn) = shape_func_c3d4(xi, eta, zeta);

                (
                    // Convert SVector to DVector using column-major slice
                    DVector::from_column_slice(n.as_slice()),
                    // IMPORTANT: Use from_column_slice because nalgebra's internal
                    // storage for SMatrix is column-major. Using from_row_slice
                    // would transpose the data and lead to a singular Jacobian.
                    DMatrix::from_column_slice(3, 4, dn.as_slice()),
                )
            }
        }
    }

    /// Computes the element stiffness matrix.
    ///
    /// # Errors
    /// Returns an error if a node is not found in the mesh or a numerical error occurs.
    pub fn compute_stiffness(
        &self,
        project: &Project,
        d_mat: &DMatrix<f64>,
        is_symmetric: bool,
    ) -> Result<DMatrix<f64>> {
        let d_static = SMatrix::<f64, 6, 6>::from_column_slice(d_mat.as_slice());

        match self {
            Element::C3D4(_, node_ids) => {
                let mut coords = SMatrix::<f64, 3, 4>::zeros();
                for (i, &node_id) in node_ids.iter().enumerate() {
                    let c = project
                        .mesh
                        .nodes
                        .get(&node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    coords.set_column(i, &nalgebra::Vector3::new(c.x, c.y, c.z));
                }

                // Fix: Pass DOF (12) explicitly as second parameter to avoid const math
                let k_static = compute_generic_stiffness::<4, 12>(
                    &d_static,
                    &coords,
                    &self.integration_points(),
                    shape_func_c3d4_static,
                    is_symmetric,
                )?;

                Ok(DMatrix::from_row_slice(12, 12, k_static.as_slice()))
            }
        }
    }

    /// Computes the element internal force vector.
    ///
    /// # Errors
    /// Returns an error if a node is not found in the mesh or numerical errors occur.
    pub fn compute_internal_force(
        &self,
        mesh: &Mesh,
        u_el: &[f64],
        d_mat: &DMatrix<f64>,
    ) -> Result<DVector<f64>> {
        match self {
            Element::C3D4(_, node_ids) => {
                let num_nodes = node_ids.len();
                let mut f_int = DVector::<f64>::zeros(num_nodes * 3);

                let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
                for (i, &node_id) in node_ids.iter().enumerate() {
                    let coords = mesh
                        .nodes
                        .get(&node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    node_coords[(0, i)] = coords.x;
                    node_coords[(1, i)] = coords.y;
                    node_coords[(2, i)] = coords.z;
                }

                for gp in self.integration_points() {
                    let (_, dn_local) =
                        self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
                    let jacobian = &dn_local * node_coords.transpose();
                    let det_j = jacobian.determinant();
                    let inv_j = jacobian
                        .try_inverse()
                        .ok_or_else(|| FerrixError::NumericalError("Singular Jacobian".into()))?;
                    let dn_global = inv_j * dn_local;
                    let weight = det_j.abs() * gp.weight;

                    // B matrix at this integration point
                    let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

                    // Stress at this integration point: sigma = D * B * u_el
                    let u_el_vec = DVector::from_column_slice(u_el);
                    let strain = &b_mat * u_el_vec;
                    let stress = d_mat * &strain;

                    // Internal force contribution: integral(B^T * sigma * dV)
                    f_int.add_assign(&(b_mat.transpose() * stress * weight));
                }

                Ok(f_int)
            }
        }
    }

    /// Calculates stress and strain at the integration point.
    ///
    /// # Errors
    /// Returns an error if a node is not found in the mesh or numerical errors occur.
    pub fn calculate_stress_strain_at_ip(
        &self,
        d_mat: &DMatrix<f64>,
        u_el: &[f64],
        mesh: &Mesh,
        ip_coords: &GaussPoint,
    ) -> Result<(DVector<f64>, DVector<f64>)> {
        let node_ids = self.get_node_ids();
        let num_nodes = node_ids.len();

        let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
        for (i, &node_id) in node_ids.iter().enumerate() {
            let coords = mesh
                .nodes
                .get(&node_id)
                .ok_or(FerrixError::NodeNotFound(node_id))?;
            node_coords[(0, i)] = coords.x;
            node_coords[(1, i)] = coords.y;
            node_coords[(2, i)] = coords.z;
        }

        let (_, dn_local) = self.shape_functions(
            ip_coords.coords[0],
            ip_coords.coords[1],
            ip_coords.coords[2],
        );

        let jacobian = &dn_local * node_coords.transpose();

        let det = jacobian.determinant();
        if det.abs() < 1e-10 {
            println!("Singular Jacobian det={det}");
            println!("Node IDs: {node_ids:?}");
            for id in node_ids {
                println!("Node {}: {:?}", id, mesh.nodes.get(id));
            }
        }

        let inv_j = jacobian
            .try_inverse()
            .ok_or_else(|| FerrixError::NumericalError("Singular Jacobian".into()))?;
        let dn_global = inv_j * dn_local;

        // Fix: Call internal build_b_matrix
        let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

        let u_el_vec = DVector::from_column_slice(u_el);
        let strain = &b_mat * u_el_vec;
        let stress = d_mat * &strain;

        Ok((strain, stress))
    }
}

fn shape_func_c3d4_static(xi: f64, eta: f64, zeta: f64) -> SMatrix<f64, 3, 4> {
    let (_, dn) = shape_func_c3d4(xi, eta, zeta);
    dn
}

// Fix: Second const param DOF to avoid { 3 * N }
/// Computes the generic stiffness matrix.
///
/// # Errors
/// Returns an error if the Jacobian is singular.
pub fn compute_generic_stiffness<const N: usize, const DOF: usize>(
    d_mat: &SMatrix<f64, 6, 6>,
    node_coords: &SMatrix<f64, 3, N>,
    integration_points: &[GaussPoint],
    shape_fn_derivatives: fn(f64, f64, f64) -> SMatrix<f64, 3, N>,
    is_symmetric: bool,
) -> Result<SMatrix<f64, DOF, DOF>> {
    let mut k_el = SMatrix::<f64, DOF, DOF>::zeros();

    for gp in integration_points {
        let dn_local = shape_fn_derivatives(gp.coords[0], gp.coords[1], gp.coords[2]);
        let jacobian = dn_local * node_coords.transpose();
        let det_j = jacobian.determinant();
        let inv_j = jacobian
            .try_inverse()
            .ok_or_else(|| FerrixError::NumericalError("Singular Jacobian".into()))?;
        let dn_global = inv_j * dn_local;
        let weight = det_j.abs() * gp.weight;

        for i in 0..N {
            let bi = build_b_block_static::<N>(&dn_global, i);
            let bit_d = bi.transpose() * d_mat;
            let start_j = if is_symmetric { i } else { 0 };

            for j in start_j..N {
                let bj = build_b_block_static::<N>(&dn_global, j);
                let k_block = (bit_d * bj) * weight;

                k_el.fixed_view_mut::<3, 3>(i * 3, j * 3)
                    .add_assign(k_block);
                if is_symmetric && i != j {
                    k_el.fixed_view_mut::<3, 3>(j * 3, i * 3)
                        .add_assign(k_block.transpose());
                }
            }
        }
    }
    Ok(k_el)
}

fn build_b_block_static<const N: usize>(
    dn_global: &SMatrix<f64, 3, N>,
    node_idx: usize,
) -> SMatrix<f64, 6, 3> {
    let mut b = SMatrix::<f64, 6, 3>::zeros();
    let dx = dn_global[(0, node_idx)];
    let dy = dn_global[(1, node_idx)];
    let dz = dn_global[(2, node_idx)];
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

// Re-implemented to fix E0425
fn build_b_matrix_internal(dn_global: &DMatrix<f64>, num_nodes: usize) -> DMatrix<f64> {
    let mut b = DMatrix::<f64>::zeros(6, num_nodes * 3);
    for i in 0..num_nodes {
        let c = i * 3;
        let dx = dn_global[(0, i)];
        let dy = dn_global[(1, i)];
        let dz = dn_global[(2, i)];
        b[(0, c)] = dx;
        b[(1, c + 1)] = dy;
        b[(2, c + 2)] = dz;
        b[(3, c)] = dy;
        b[(3, c + 1)] = dx;
        b[(4, c + 1)] = dz;
        b[(4, c + 2)] = dy;
        b[(5, c)] = dz;
        b[(5, c + 2)] = dx;
    }
    b
}

#[derive(Debug, Clone, Copy)]
pub struct GaussPoint {
    pub coords: [f64; 3],
    pub weight: f64,
}
