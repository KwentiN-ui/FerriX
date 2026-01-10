use crate::solver::assembler::build_b_matrix;
use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::mesh::Mesh;
use nalgebra::{DMatrix, DVector};
use std::error::Error;
use std::str::FromStr;

use strum_macros::{EnumDiscriminants, EnumString};

/// Supported element types (<https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node194.html>).
/// `strum` automatically generates a String-enum `ElementType` based on these definitions.
#[derive(EnumDiscriminants, Debug, Clone)]
#[strum_discriminants(derive(Hash, EnumString))]
#[strum_discriminants(name(ElementType))]
pub enum Element {
    // General 3D-Solids
    /// 4-node linear tetrahedral element
    C3D4(ElementId, [NodeId; 4]),
    // TODO: 3D 20-node quadratic isoparametric element
    // C3D20(ElementId, [NodeId; 20]),
}

impl Element {
    pub fn parse_type_str_from_line(line: &str) -> Result<String, Box<dyn Error>> {
        Ok(line
            .split(',')
            .map(str::trim)
            .nth(1)
            .ok_or("Invalid element definition on line {line_nr}")?
            .split('=')
            .next_back()
            .ok_or("Invalid element definition on line {line_nr}")?
            .to_string())
    }
    /// Create an Element from a line. This function panics if it's not able to do so.
    pub fn parse_line(type_name: &str, line: &str) -> Self {
        let nums: Vec<usize> = line
            .split(',')
            .map(|s| s.trim().parse().expect("Integer conversion failed"))
            .collect();

        let (&id, nodes_usize) = nums.split_first().expect("Line empty");
        let id = ElementId(id);
        let nodes: Vec<NodeId> = nodes_usize.iter().map(|&n| NodeId(n)).collect();

        // String -> ElementType
        let elem_type = ElementType::from_str(type_name)
            .unwrap_or_else(|_| panic!("Unknown element definition: {type_name}"));

        // Local Macro for array casting
        macro_rules! to_arr {
            ($n:expr) => {
                nodes
                    .try_into()
                    .expect(concat!("Wrong node count for ", stringify!($n)))
            };
        }

        match elem_type {
            ElementType::C3D4 => Element::C3D4(id, to_arr!(C3D4)),
        }
    }

    pub fn get_id(&self) -> ElementId {
        match self {
            Element::C3D4(id, _) => *id,
        }
    }

    /// Get global Node IDs
    pub fn get_node_ids(&self) -> &[NodeId] {
        match self {
            Element::C3D4(_, n) => n,
        }
    }

    pub fn integration_points(&self) -> Vec<GaussPoint> {
        match self {
            Element::C3D4(..) => c3d4_gauss(),
        }
    }

    /// Retrieves the element types shape functions, and shape function derivatives.
    pub fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        match self {
            Element::C3D4(..) => shape_func_c3d4(xi, eta, zeta),
        }
    }

    pub fn calculate_stress_strain_at_ip(
        &self,
        d_mat: &DMatrix<f64>,
        u_el: &[f64], // Nodal displacements for this element
        mesh: &Mesh,
        ip_coords: &GaussPoint,
    ) -> (DVector<f64>, DVector<f64>) {
        let node_ids = self.get_node_ids();
        let num_nodes = node_ids.len();

        let mut node_coords = DMatrix::<f64>::zeros(3, num_nodes);
        for (i, &node_id) in node_ids.iter().enumerate() {
            let coords = mesh.nodes.get(&node_id).expect("Node not found");
            node_coords[(0, i)] = coords.x;
            node_coords[(1, i)] = coords.y;
            node_coords[(2, i)] = coords.z;
        }

        // 2. Compute mapping and B-Matrix
        let (_, dn_local) = self.shape_functions(
            ip_coords.coords[0],
            ip_coords.coords[1],
            ip_coords.coords[2],
        );

        // Jacobian: J = dN * X^T (3x3 matrix)
        let jacobian = &dn_local * node_coords.transpose();
        let inv_j = jacobian.try_inverse().expect("Singular Jacobian");
        let dn_global = inv_j * dn_local;
        let b_mat = build_b_matrix(&dn_global, num_nodes);

        // 3. Replace u_el_array with DVector
        let u_el_vec = DVector::from_column_slice(u_el);

        // 4. Calculate strain and stress using standard matrix-vector multiplication
        // epsilon = B * u
        let strain = &b_mat * u_el_vec;
        // sigma = D * epsilon
        let stress = (d_mat * &strain)
            // needed for type safe declaration of column-vector
            .column(0)
            .into();

        (strain, stress)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GaussPoint {
    pub coords: [f64; 3], // xi, eta, zeta
    pub weight: f64,
}

use crate::solver::mesh_lib::elements::c3d4::{c3d4_gauss, shape_func_c3d4};
