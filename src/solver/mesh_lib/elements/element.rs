//! Base element definitions and generic element logic.
//!
//! This module defines the `Element` enum, which encapsulates all supported
//! finite element types, and providing methods for stiffness and force calculations.

use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::{ElementId, NodeId};
use crate::solver::material::Material;
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::project::Project;
use nalgebra::{DMatrix, DVector, SMatrix};
use std::ops::AddAssign;
use std::str::FromStr;

use crate::solver::mesh_lib::elements::c3d4::{c3d4_gauss, shape_func_c3d4};
use strum_macros::{EnumDiscriminants, EnumString};

/// Supported finite element types.
#[derive(EnumDiscriminants, Debug, Clone)]
#[strum_discriminants(derive(Hash, EnumString))]
#[strum_discriminants(name(ElementType))]
pub enum Element {
    /// A 4-node linear tetrahedron (First-order 3D element).
    C3D4(ElementId, [NodeId; 4]),
}

impl Element {
    /// Parses the element type name from a keyword line (e.g., `*ELEMENT, TYPE=C3D4`).
    ///
    /// # Errors
    /// Returns `FerrixError::ParseError` if the line is malformed.
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

    /// Parses an element's ID and connectivity from a data line.
    ///
    /// # Errors
    /// Returns an error if the input line is malformed or the element type is unknown.
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

    /// Returns the unique ID of the element.
    #[must_use]
    pub fn get_id(&self) -> ElementId {
        match self {
            Element::C3D4(id, _) => *id,
        }
    }

    /// Returns a slice of node IDs that form the element's connectivity.
    #[must_use]
    pub fn get_node_ids(&self) -> &[NodeId] {
        match self {
            Element::C3D4(_, n) => n,
        }
    }

    /// Returns the Gauss integration points for this element type.
    #[must_use]
    pub fn integration_points(&self) -> Vec<GaussPoint> {
        match self {
            Element::C3D4(..) => c3d4_gauss(),
        }
    }

    /// Computes shape functions (N) and their derivatives (dN) at local coordinates (xi, eta, zeta).
    #[must_use]
    pub fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        match self {
            Element::C3D4(..) => {
                let (n, dn) = shape_func_c3d4(xi, eta, zeta);
                (
                    DVector::from_column_slice(n.as_slice()),
                    DMatrix::from_column_slice(3, 4, dn.as_slice()),
                )
            }
        }
    }

    /// Computes the element's local stiffness matrix (`k_el`).
    ///
    /// # Errors
    /// Returns an error if a node is missing from the mesh or if the Jacobian is singular.
    pub fn compute_stiffness(
        &self,
        project: &Project,
        d_mat: &DMatrix<f64>,
        is_symmetric: bool,
        u_el: Option<&[f64]>,
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

                    let mut pos = nalgebra::Vector3::new(c.x, c.y, c.z);
                    if let Some(u) = u_el {
                        pos[0] += u[i * 3];
                        pos[1] += u[i * 3 + 1];
                        pos[2] += u[i * 3 + 2];
                    }
                    coords.set_column(i, &pos);
                }

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

    /// Computes the element's internal force vector.
    ///
    /// # Errors
    /// Returns an error if node coordinates are missing or if numerical issues arise during integration.
    pub fn compute_internal_force(
        &self,
        mesh: &Mesh,
        material: &dyn Material,
        u_el: &[f64],
        node_temps: &[f64],
        u_conf: Option<&[f64]>,
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

                    let mut pos = nalgebra::Vector3::new(coords.x, coords.y, coords.z);
                    if let Some(u) = u_conf {
                        pos[0] += u[i * 3];
                        pos[1] += u[i * 3 + 1];
                        pos[2] += u[i * 3 + 2];
                    }
                    node_coords.set_column(i, &pos);
                }

                for gp in self.integration_points() {
                    let (n_local, dn_local) =
                        self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
                    let jacobian = &dn_local * node_coords.transpose();
                    let det_j = jacobian.determinant();
                    let inv_j = jacobian
                        .try_inverse()
                        .ok_or_else(|| FerrixError::NumericalError("Singular Jacobian".into()))?;
                    let dn_global = inv_j * dn_local;
                    let weight = det_j.abs() * gp.weight;

                    let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

                    // Interpolate temperature at integration point
                    let mut t_ip = 0.0;
                    for i in 0..num_nodes {
                        t_ip += n_local[i] * node_temps[i];
                    }

                    let d_mat = material.build_elastic_d_matrix(t_ip)?;

                    let u_el_vec = DVector::from_column_slice(u_el);
                    let mut strain = &b_mat * u_el_vec;

                    // Subtract thermal strain if expansion is defined
                    if let Some(alpha) = material.thermal_expansion(t_ip) {
                        let t_ref = material.reference_temperature();
                        let delta_t = t_ip - t_ref;
                        strain[0] -= alpha * delta_t;
                        strain[1] -= alpha * delta_t;
                        strain[2] -= alpha * delta_t;
                    }

                    let stress = d_mat * &strain;

                    f_int.add_assign(&(b_mat.transpose() * stress * weight));
                }

                Ok(f_int)
            }
        }
    }

    /// Computes the element's thermal force vector.
    ///
    /// # Errors
    /// Returns an error if node coordinates or temperatures are missing.
    pub fn compute_thermal_force(
        &self,
        project: &Project,
        material: &dyn Material,
        node_temps: &[f64],
    ) -> Result<DVector<f64>> {
        match self {
            Element::C3D4(_, node_ids) => {
                let num_nodes = node_ids.len();
                let mut f_th = DVector::<f64>::zeros(num_nodes * 3);

                let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
                for (i, &node_id) in node_ids.iter().enumerate() {
                    let coords = project
                        .mesh
                        .nodes
                        .get(&node_id)
                        .ok_or(FerrixError::NodeNotFound(node_id))?;
                    node_coords
                        .set_column(i, &nalgebra::Vector3::new(coords.x, coords.y, coords.z));
                }

                for gp in self.integration_points() {
                    let (n_local, dn_local) =
                        self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);

                    // Interpolate temperature at integration point
                    let mut t_ip = 0.0;
                    for i in 0..num_nodes {
                        t_ip += n_local[i] * node_temps[i];
                    }

                    let Some(alpha) = material.thermal_expansion(t_ip) else {
                        continue;
                    };
                    let t_ref = material.reference_temperature();
                    let d_mat = material.build_elastic_d_matrix(t_ip)?;

                    let jacobian = &dn_local * node_coords.transpose();
                    let det_j = jacobian.determinant();
                    let inv_j = jacobian
                        .try_inverse()
                        .ok_or_else(|| FerrixError::NumericalError("Singular Jacobian".into()))?;
                    let dn_global = inv_j * dn_local;
                    let weight = det_j.abs() * gp.weight;

                    let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

                    // Thermal strain vector: [alpha*deltaT, alpha*deltaT, alpha*deltaT, 0, 0, 0]
                    let delta_t = t_ip - t_ref;
                    let mut strain_th = DVector::<f64>::zeros(6);
                    strain_th[0] = alpha * delta_t;
                    strain_th[1] = alpha * delta_t;
                    strain_th[2] = alpha * delta_t;

                    let stress_th = d_mat * strain_th;
                    f_th.add_assign(&(b_mat.transpose() * stress_th * weight));
                }

                Ok(f_th)
            }
        }
    }

    /// Calculates stress and strain vectors at a specific integration point.
    ///
    /// # Errors
    /// Returns an error if numerical errors occur during Jacobian inversion.
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
        let inv_j = jacobian
            .try_inverse()
            .ok_or_else(|| FerrixError::NumericalError("Singular Jacobian".into()))?;
        let dn_global = inv_j * dn_local;

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

/// Generic routine for computing the stiffness matrix for any element type.
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

/// Represents a single integration point (Gauss point) within an element.
#[derive(Debug, Clone, Copy)]
pub struct GaussPoint {
    /// Local (isoparametric) coordinates of the point.
    pub coords: [f64; 3],
    /// Integration weight associated with the point.
    pub weight: f64,
}
