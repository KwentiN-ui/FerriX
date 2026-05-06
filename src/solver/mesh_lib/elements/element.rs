//! Base element definitions and generic element logic.
//!
//! This module defines the `Element` enum and `FiniteElement` trait, which encapsulates all supported
//! finite element types, and providing methods for stiffness and force calculations.

use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::{ElementId, NodeId};
use crate::solver::material::{ElementMaterialState, Material, MaterialPointState};
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::project::Project;
use nalgebra::{DMatrix, DVector};
use std::ops::AddAssign;
use std::str::FromStr;

use crate::solver::mesh_lib::elements::c3d4::C3D4;
use crate::solver::mesh_lib::elements::c3d6::C3D6;
use crate::solver::mesh_lib::elements::c3d10::C3D10;
use crate::solver::mesh_lib::elements::c3d20::C3D20;
use strum_macros::{EnumDiscriminants, EnumString};

/// Trait defining the mathematical and topological properties of a finite element.
pub trait FiniteElement: std::fmt::Debug + Send + Sync {
    /// Returns the unique ID of the element.
    fn id(&self) -> ElementId;
    /// Returns a slice of node IDs that form the element's connectivity.
    fn nodes(&self) -> &[NodeId];
    /// Returns the number of nodes in this element.
    fn num_nodes(&self) -> usize;
    /// Returns the VTK cell type code for this element.
    fn vtk_cell_type(&self) -> u8;
    /// Returns the Gauss integration points for this element type.
    fn integration_points(&self) -> &'static [GaussPoint];
    /// Computes shape functions (N) and their derivatives (dN) at local coordinates (xi, eta, zeta).
    fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>);
    /// Returns the local coordinates of the nodes.
    // fn node_local_coords(&self) -> Vec<[f64; 3]>;
    fn node_local_coords(&self) -> &'static [[f64; 3]];
}

/// Supported finite element types.
#[derive(EnumDiscriminants, Debug, Clone)]
#[strum_discriminants(derive(Hash, EnumString))]
#[strum_discriminants(name(ElementType))]
pub enum Element {
    /// A 4-node linear tetrahedron (First-order 3D element).
    C3D4(C3D4),
    /// A 6-node linear wedge (First-order 3D element).
    C3D6(C3D6),
    /// A 10-node quadratic tetrahedron (Second-order 3D element).
    C3D10(C3D10),
    /// A 20-node quadratic brick (Second-order 3D element).
    C3D20(C3D20),
}

impl Element {
    /// Returns a reference to the inner `FiniteElement` implementation.
    #[must_use]
    pub fn inner(&self) -> &dyn FiniteElement {
        match self {
            Element::C3D4(e) => e,
            Element::C3D6(e) => e,
            Element::C3D10(e) => e,
            Element::C3D20(e) => e,
        }
    }

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
                Ok(Element::C3D4(C3D4 {
                    id,
                    nodes: nodes_arr,
                }))
            }
            ElementType::C3D6 => {
                let nodes_arr: [NodeId; 6] =
                    nodes.try_into().map_err(|_| FerrixError::ParseError {
                        line: 0,
                        message: "Wrong node count for C3D6".into(),
                    })?;
                Ok(Element::C3D6(C3D6 {
                    id,
                    nodes: nodes_arr,
                }))
            }
            ElementType::C3D10 => {
                let nodes_arr: [NodeId; 10] =
                    nodes.try_into().map_err(|_| FerrixError::ParseError {
                        line: 0,
                        message: "Wrong node count for C3D10".into(),
                    })?;
                Ok(Element::C3D10(C3D10 {
                    id,
                    nodes: nodes_arr,
                }))
            }
            ElementType::C3D20 => {
                let nodes_arr: [NodeId; 20] =
                    nodes.try_into().map_err(|_| FerrixError::ParseError {
                        line: 0,
                        message: "Wrong node count for C3D20".into(),
                    })?;
                Ok(Element::C3D20(C3D20 {
                    id,
                    nodes: nodes_arr,
                }))
            }
        }
    }

    /// Returns the unique ID of the element.
    #[must_use]
    pub fn get_id(&self) -> ElementId {
        self.inner().id()
    }

    /// Returns a slice of node IDs that form the element's connectivity.
    #[must_use]
    pub fn get_node_ids(&self) -> &[NodeId] {
        self.inner().nodes()
    }

    /// Returns the Gauss integration points for this element type.
    #[must_use]
    pub fn integration_points(&self) -> &'static [GaussPoint] {
        self.inner().integration_points()
    }

    /// Computes shape functions (N) and their derivatives (dN) at local coordinates (xi, eta, zeta).
    #[must_use]
    pub fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        self.inner().shape_functions(xi, eta, zeta)
    }

    /// Returns the local coordinates of the nodes.
    #[must_use]
    pub fn node_local_coords(&self) -> &'static [[f64; 3]] {
        self.inner().node_local_coords()
    }

    /// Computes the element's local stiffness matrix (`k_el`).
    ///
    /// # Errors
    /// Returns an error if a node is missing from the mesh or if the Jacobian is singular.
    pub fn compute_stiffness(
        &self,
        project: &Project,
        d_mat: &DMatrix<f64>,
        u_el: Option<&[f64]>,
    ) -> Result<DMatrix<f64>> {
        let node_ids = self.get_node_ids();
        let num_nodes = node_ids.len();

        let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
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
            node_coords.set_column(i, &pos);
        }

        let mut k_el = DMatrix::<f64>::zeros(num_nodes * 3, num_nodes * 3);

        for gp in self.integration_points() {
            let (_, dn_local) = self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
            let jacobian = &dn_local * node_coords.transpose();
            let det_j = jacobian.determinant();
            let inv_j = jacobian.try_inverse().ok_or_else(|| {
                FerrixError::NumericalError(format!(
                    "Singular Jacobian in element {}",
                    self.get_id().0
                ))
            })?;
            let dn_global = inv_j * dn_local;
            let weight = det_j.abs() * gp.weight;

            let b_mat = build_b_matrix_internal(&dn_global, num_nodes);
            k_el.add_assign(&(b_mat.transpose() * d_mat * b_mat * weight));
        }

        Ok(k_el)
    }

    /// Computes element stiffness and updated SDVs.
    ///
    /// # Errors
    /// Returns `FerrixError` if numerical errors occur or material update fails.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn compute_stiffness_sdv(
        &self,
        project: &Project,
        material: &dyn Material,
        u_el: Option<&[f64]>,
        node_temps_initial: Option<&[f64]>,
        node_temps_current: Option<&[f64]>,
        material_states_old: Option<&ElementMaterialState>,
        dtime: f64,
    ) -> Result<(DMatrix<f64>, Option<ElementMaterialState>)> {
        let node_ids = self.get_node_ids();
        let num_nodes = node_ids.len();
        let num_dofs = num_nodes * 3;
        let mut k_el = DMatrix::<f64>::zeros(num_dofs, num_dofs);
        let mut updated_states = Vec::new();

        let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
        for (i, &node_id) in node_ids.iter().enumerate() {
            let coords = project
                .mesh
                .nodes
                .get(&node_id)
                .ok_or(FerrixError::NodeNotFound(node_id))?;

            let mut pos = nalgebra::Vector3::new(coords.x, coords.y, coords.z);
            if let Some(u) = u_el {
                pos[0] += u[i * 3];
                pos[1] += u[i * 3 + 1];
                pos[2] += u[i * 3 + 2];
            }
            node_coords.set_column(i, &pos);
        }

        let empty_state = MaterialPointState::default();

        for (ip_idx, gp) in self.integration_points().iter().enumerate() {
            let (n_local, dn_local) =
                self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
            let jacobian = &dn_local * node_coords.transpose();
            let det_j = jacobian.determinant();
            let inv_j = jacobian.try_inverse().ok_or_else(|| {
                FerrixError::NumericalError(format!(
                    "Singular Jacobian in element {}",
                    self.get_id().0
                ))
            })?;
            let dn_global = inv_j * dn_local;
            let weight = det_j.abs() * gp.weight;

            let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

            let t_curr = node_temps_current.map_or(0.0, |temps| {
                let mut t = 0.0;
                for i in 0..num_nodes {
                    t += n_local[i] * temps[i];
                }
                t
            });

            let t_init = node_temps_initial.map_or(t_curr, |temps| {
                let mut t = 0.0;
                for i in 0..num_nodes {
                    t += n_local[i] * temps[i];
                }
                t
            });

            let u_el_vec = u_el.map_or(DVector::zeros(num_dofs), DVector::from_column_slice);
            let mut strain = &b_mat * u_el_vec;

            // Subtract thermal strain delta: epsilon_th(T_curr) - epsilon_th(T_init)
            if let Some(alpha_curr) = material.thermal_expansion(t_curr) {
                let t_ref = material.reference_temperature();
                let th_strain_curr = alpha_curr * (t_curr - t_ref);

                let alpha_init = material.thermal_expansion(t_init).unwrap_or(alpha_curr);
                let th_strain_init = alpha_init * (t_init - t_ref);

                let delta_th_strain = th_strain_curr - th_strain_init;

                strain[0] -= delta_th_strain;
                strain[1] -= delta_th_strain;
                strain[2] -= delta_th_strain;
            }

            let state_old = material_states_old.map_or(&empty_state, |m| &m[ip_idx]);
            let (d_tangent, _stress, state_new) =
                material.update_state(t_curr, &strain, state_old, dtime)?;
            updated_states.push(state_new);

            k_el.add_assign(&(b_mat.transpose() * d_tangent * b_mat * weight));
        }

        let has_sdvs = material.num_state_variables() > 0;
        Ok((
            k_el,
            if has_sdvs {
                Some(ElementMaterialState {
                    ip_states: updated_states,
                })
            } else {
                None
            },
        ))
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
        node_temps_initial: &[f64],
        node_temps_current: &[f64],
        u_conf: Option<&[f64]>,
    ) -> Result<DVector<f64>> {
        // Just call the SDV version with None for old states
        let (f_int, _) = self.compute_internal_force_sdv(
            mesh,
            material,
            u_el,
            node_temps_initial,
            node_temps_current,
            None,
            0.0,
            u_conf,
        )?;
        Ok(f_int)
    }

    /// Computes internal force and updated SDVs.
    ///
    /// # Errors
    /// Returns `FerrixError` if numerical errors occur or material update fails.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn compute_internal_force_sdv(
        &self,
        mesh: &Mesh,
        material: &dyn Material,
        u_el: &[f64],
        node_temps_initial: &[f64],
        node_temps_current: &[f64],
        material_states_old: Option<&ElementMaterialState>,
        dtime: f64,
        u_conf: Option<&[f64]>,
    ) -> Result<(DVector<f64>, Option<ElementMaterialState>)> {
        let node_ids = self.get_node_ids();
        let num_nodes = node_ids.len();
        let num_dofs = num_nodes * 3;
        let mut f_int = DVector::<f64>::zeros(num_dofs);
        let mut updated_states = Vec::new();

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

        let empty_state = MaterialPointState::default();

        for (ip_idx, gp) in self.integration_points().iter().enumerate() {
            let (n_local, dn_local) =
                self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
            let jacobian = &dn_local * node_coords.transpose();
            let det_j = jacobian.determinant();
            let inv_j = jacobian.try_inverse().ok_or_else(|| {
                FerrixError::NumericalError(format!(
                    "Singular Jacobian in element {}",
                    self.get_id().0
                ))
            })?;
            let dn_global = inv_j * dn_local;
            let weight = det_j.abs() * gp.weight;

            let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

            // Interpolate temperatures at integration point
            let mut t_curr = 0.0;
            let mut t_init = 0.0;
            for i in 0..num_nodes {
                t_curr += n_local[i] * node_temps_current[i];
                t_init += n_local[i] * node_temps_initial[i];
            }

            let _d_mat_check = material.build_elastic_d_matrix(t_curr)?;

            let u_el_vec = DVector::from_column_slice(u_el);
            let mut strain = &b_mat * u_el_vec;

            // Subtract thermal strain delta: epsilon_th(T_curr) - epsilon_th(T_init)
            if let Some(alpha_curr) = material.thermal_expansion(t_curr) {
                let t_ref = material.reference_temperature();
                let th_strain_curr = alpha_curr * (t_curr - t_ref);

                let alpha_init = material.thermal_expansion(t_init).unwrap_or(alpha_curr);
                let th_strain_init = alpha_init * (t_init - t_ref);

                let delta_th_strain = th_strain_curr - th_strain_init;

                strain[0] -= delta_th_strain;
                strain[1] -= delta_th_strain;
                strain[2] -= delta_th_strain;
            }

            let state_old = material_states_old.map_or(&empty_state, |m| &m[ip_idx]);
            let (_d_tangent, stress, state_new) =
                material.update_state(t_curr, &strain, state_old, dtime)?;
            updated_states.push(state_new);

            f_int.add_assign(&(b_mat.transpose() * stress * weight));
        }

        let has_sdvs = material.num_state_variables() > 0;
        Ok((
            f_int,
            if has_sdvs {
                Some(ElementMaterialState {
                    ip_states: updated_states,
                })
            } else {
                None
            },
        ))
    }

    /// Computes the element's thermal force vector.
    ///
    /// # Errors
    /// Returns an error if node coordinates or temperatures are missing.
    pub fn compute_thermal_force(
        &self,
        project: &Project,
        material: &dyn Material,
        node_temps_initial: &[f64],
        node_temps_current: &[f64],
    ) -> Result<DVector<f64>> {
        let node_ids = self.get_node_ids();
        let num_nodes = node_ids.len();
        let mut f_th = DVector::<f64>::zeros(num_nodes * 3);

        let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
        for (i, &node_id) in node_ids.iter().enumerate() {
            let coords = project
                .mesh
                .nodes
                .get(&node_id)
                .ok_or(FerrixError::NodeNotFound(node_id))?;
            node_coords.set_column(i, &nalgebra::Vector3::new(coords.x, coords.y, coords.z));
        }

        for gp in self.integration_points() {
            let (n_local, dn_local) =
                self.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);

            // Interpolate temperatures at integration point
            let mut t_curr = 0.0;
            let mut t_init = 0.0;
            for i in 0..num_nodes {
                t_curr += n_local[i] * node_temps_current[i];
                t_init += n_local[i] * node_temps_initial[i];
            }

            let Some(alpha_curr) = material.thermal_expansion(t_curr) else {
                continue;
            };
            let t_ref = material.reference_temperature();
            let th_strain_curr = alpha_curr * (t_curr - t_ref);

            let alpha_init = material.thermal_expansion(t_init).unwrap_or(alpha_curr);
            let th_strain_init = alpha_init * (t_init - t_ref);

            let delta_th_strain = th_strain_curr - th_strain_init;

            let d_mat = material.build_elastic_d_matrix(t_curr)?;

            let jacobian = &dn_local * node_coords.transpose();
            let det_j = jacobian.determinant();
            let inv_j = jacobian.try_inverse().ok_or_else(|| {
                FerrixError::NumericalError(format!(
                    "Singular Jacobian in element {}",
                    self.get_id().0
                ))
            })?;
            let dn_global = inv_j * dn_local;
            let weight = det_j.abs() * gp.weight;

            let b_mat = build_b_matrix_internal(&dn_global, num_nodes);

            // Thermal strain delta vector
            let mut strain_th = DVector::<f64>::zeros(6);
            strain_th[0] = delta_th_strain;
            strain_th[1] = delta_th_strain;
            strain_th[2] = delta_th_strain;

            let stress_th = d_mat * strain_th;
            f_th.add_assign(&(b_mat.transpose() * stress_th * weight));
        }

        Ok(f_th)
    }

    /// Calculates stress and strain vectors at a specific local coordinates.
    ///
    /// # Errors
    /// Returns an error if numerical errors occur during Jacobian inversion.
    pub fn calculate_stress_strain_at_local_coords(
        &self,
        d_mat: &DMatrix<f64>,
        u_el: &[f64],
        mesh: &Mesh,
        xi: f64,
        eta: f64,
        zeta: f64,
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

        let (_, dn_local) = self.shape_functions(xi, eta, zeta);

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
        self.calculate_stress_strain_at_local_coords(
            d_mat,
            u_el,
            mesh,
            ip_coords.coords[0],
            ip_coords.coords[1],
            ip_coords.coords[2],
        )
    }
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
