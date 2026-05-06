//! C3D6 element implementation (Linear Wedge).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector};

/// Constant data for C3D6 integration points.
static C3D6_N: [[f64; 6]; 2] = [
    c3d6_math(
        0.333_333_333_333_333,
        0.333_333_333_333_333,
        -0.577_350_269_189_626,
    )
    .0,
    c3d6_math(
        0.333_333_333_333_333,
        0.333_333_333_333_333,
        0.577_350_269_189_626,
    )
    .0,
];
static C3D6_DN: [[f64; 18]; 2] = [
    c3d6_math(
        0.333_333_333_333_333,
        0.333_333_333_333_333,
        -0.577_350_269_189_626,
    )
    .1,
    c3d6_math(
        0.333_333_333_333_333,
        0.333_333_333_333_333,
        0.577_350_269_189_626,
    )
    .1,
];

/// Gaussian integration points for C3D6.
const C3D6_GAUSS: [GaussPoint; 2] = [
    GaussPoint {
        coords: [
            0.333_333_333_333_333,
            0.333_333_333_333_333,
            -0.577_350_269_189_626,
        ],
        weight: 0.5,
        n: &C3D6_N[0],
        dn: &C3D6_DN[0],
    },
    GaussPoint {
        coords: [
            0.333_333_333_333_333,
            0.333_333_333_333_333,
            0.577_350_269_189_626,
        ],
        weight: 0.5,
        n: &C3D6_N[1],
        dn: &C3D6_DN[1],
    },
];

const C3D6_LOCAL_COORDS: [[f64; 3]; 6] = [
    [0.0, 0.0, -1.0], // Node 1
    [1.0, 0.0, -1.0], // Node 2
    [0.0, 1.0, -1.0], // Node 3
    [0.0, 0.0, 1.0],  // Node 4
    [1.0, 0.0, 1.0],  // Node 5
    [0.0, 1.0, 1.0],  // Node 6
];

/// Linear wedge (C3D6).
#[derive(Debug, Clone)]
pub struct C3D6 {
    pub id: ElementId,
    pub nodes: [NodeId; 6],
}

impl FiniteElement for C3D6 {
    fn id(&self) -> ElementId {
        self.id
    }

    fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    fn num_nodes(&self) -> usize {
        6
    }

    fn vtk_cell_type(&self) -> u8 {
        13 // VTK_WEDGE
    }

    fn integration_points(&self) -> &'static [GaussPoint] {
        &C3D6_GAUSS
    }

    fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        let (n, dn) = c3d6_math(xi, eta, zeta);
        (
            DVector::from_column_slice(&n),
            DMatrix::from_column_slice(3, 6, &dn),
        )
    }

    fn node_local_coords(&self) -> &'static [[f64; 3]] {
        &C3D6_LOCAL_COORDS
    }
}

/// Mathematical definition of C3D6 shape functions and their derivatives.
///
/// Returns (N, dN) where N is a 6-element array and dN is a flattened 3x6 matrix
/// in COLUMN-MAJOR order.
const fn c3d6_math(xi: f64, et: f64, ze: f64) -> ([f64; 6], [f64; 18]) {
    let a = 1.0 - xi - et;

    let n = [
        0.5 * a * (1.0 - ze),
        0.5 * xi * (1.0 - ze),
        0.5 * et * (1.0 - ze),
        0.5 * a * (1.0 + ze),
        0.5 * xi * (1.0 + ze),
        0.5 * et * (1.0 + ze),
    ];

    let dn = [
        // Node 1
        -0.5 * (1.0 - ze), // dN/dxi
        -0.5 * (1.0 - ze), // dN/deta
        -0.5 * a,          // dN/dzeta
        // Node 2
        0.5 * (1.0 - ze),
        0.0,
        -0.5 * xi,
        // Node 3
        0.0,
        0.5 * (1.0 - ze),
        -0.5 * et,
        // Node 4
        -0.5 * (1.0 + ze),
        -0.5 * (1.0 + ze),
        0.5 * a,
        // Node 5
        0.5 * (1.0 + ze),
        0.0,
        0.5 * xi,
        // Node 6
        0.0,
        0.5 * (1.0 + ze),
        0.5 * et,
    ];

    (n, dn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_func_sum() {
        let (n, _) = c3d6_math(0.25, 0.25, 0.0);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        let (n, _) = c3d6_math(0.1, 0.2, 0.5);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shape_func_nodes() {
        for (i, coord) in C3D6_LOCAL_COORDS.iter().enumerate() {
            let (n, _) = c3d6_math(coord[0], coord[1], coord[2]);
            assert!(
                (n[i] - 1.0).abs() < 1e-12,
                "Node {} failed: n[{}] = {}",
                i + 1,
                i,
                n[i]
            );
            for (j, nj) in n.iter().enumerate() {
                if i != j {
                    assert!(
                        nj.abs() < 1e-12,
                        "Node {} failed at n[{}]: n[{}] = {}",
                        i + 1,
                        j,
                        j,
                        nj
                    );
                }
            }
        }
    }
}
