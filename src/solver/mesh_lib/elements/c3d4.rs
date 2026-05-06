//! C3D4 element implementation (Linear Tetrahedron).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector};

/// Constant data for C3D4 integration points.
static C3D4_N: [[f64; 4]; 1] = [c3d4_math(0.25, 0.25, 0.25).0];
static C3D4_DN: [[f64; 12]; 1] = [c3d4_math(0.25, 0.25, 0.25).1];

const C3D4_GAUSS: [GaussPoint; 1] = [GaussPoint {
    coords: [0.25, 0.25, 0.25],
    weight: 1.0 / 6.0,
    n: &C3D4_N[0],
    dn: &C3D4_DN[0],
}];

const C3D4_LOCAL_COORDS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// Linear tetrahedron (C3D4).
#[derive(Debug, Clone)]
pub struct C3D4 {
    pub id: ElementId,
    pub nodes: [NodeId; 4],
}

impl FiniteElement for C3D4 {
    fn id(&self) -> ElementId {
        self.id
    }

    fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    fn num_nodes(&self) -> usize {
        4
    }

    fn vtk_cell_type(&self) -> u8 {
        10 // VTK_TETRA
    }

    fn integration_points(&self) -> &'static [GaussPoint] {
        &C3D4_GAUSS
    }

    fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        let (n, dn) = c3d4_math(xi, eta, zeta);
        (
            DVector::from_column_slice(&n),
            DMatrix::from_column_slice(3, 4, &dn),
        )
    }

    fn node_local_coords(&self) -> &'static [[f64; 3]] {
        &C3D4_LOCAL_COORDS
    }
}

/// Mathematical definition of C3D4 shape functions and their derivatives.
///
/// Returns (N, dN) where N is a 4-element array and dN is a flattened 3x4 matrix
/// in COLUMN-MAJOR order.
const fn c3d4_math(xi: f64, eta: f64, zeta: f64) -> ([f64; 4], [f64; 12]) {
    let n = [1.0 - xi - eta - zeta, xi, eta, zeta];

    // Column-major order: [dN1/dxi, dN1/deta, dN1/dzeta, dN2/dxi, ...]
    let dn = [
        -1.0, -1.0, -1.0, // Node 1
        1.0, 0.0, 0.0, // Node 2
        0.0, 1.0, 0.0, // Node 3
        0.0, 0.0, 1.0, // Node 4
    ];
    (n, dn)
}
