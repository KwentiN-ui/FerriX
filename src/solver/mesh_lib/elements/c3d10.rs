//! C3D10 element implementation (Quadratic Tetrahedron).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector};

/// Gaussian integration points for C3D10.
const C3D10_GAUSS_COORDS: [[f64; 3]; 4] = [
    [
        0.138_196_601_125_011,
        0.138_196_601_125_011,
        0.138_196_601_125_011,
    ],
    [
        0.585_410_196_624_968,
        0.138_196_601_125_011,
        0.138_196_601_125_011,
    ],
    [
        0.138_196_601_125_011,
        0.585_410_196_624_968,
        0.138_196_601_125_011,
    ],
    [
        0.138_196_601_125_011,
        0.138_196_601_125_011,
        0.585_410_196_624_968,
    ],
];

/// Constant data for C3D10 integration points.
static C3D10_N: [[f64; 10]; 4] = [
    c3d10_math(
        C3D10_GAUSS_COORDS[0][0],
        C3D10_GAUSS_COORDS[0][1],
        C3D10_GAUSS_COORDS[0][2],
    )
    .0,
    c3d10_math(
        C3D10_GAUSS_COORDS[1][0],
        C3D10_GAUSS_COORDS[1][1],
        C3D10_GAUSS_COORDS[1][2],
    )
    .0,
    c3d10_math(
        C3D10_GAUSS_COORDS[2][0],
        C3D10_GAUSS_COORDS[2][1],
        C3D10_GAUSS_COORDS[2][2],
    )
    .0,
    c3d10_math(
        C3D10_GAUSS_COORDS[3][0],
        C3D10_GAUSS_COORDS[3][1],
        C3D10_GAUSS_COORDS[3][2],
    )
    .0,
];
static C3D10_DN: [[f64; 30]; 4] = [
    c3d10_math(
        C3D10_GAUSS_COORDS[0][0],
        C3D10_GAUSS_COORDS[0][1],
        C3D10_GAUSS_COORDS[0][2],
    )
    .1,
    c3d10_math(
        C3D10_GAUSS_COORDS[1][0],
        C3D10_GAUSS_COORDS[1][1],
        C3D10_GAUSS_COORDS[1][2],
    )
    .1,
    c3d10_math(
        C3D10_GAUSS_COORDS[2][0],
        C3D10_GAUSS_COORDS[2][1],
        C3D10_GAUSS_COORDS[2][2],
    )
    .1,
    c3d10_math(
        C3D10_GAUSS_COORDS[3][0],
        C3D10_GAUSS_COORDS[3][1],
        C3D10_GAUSS_COORDS[3][2],
    )
    .1,
];

const C3D10_GAUSS: [GaussPoint; 4] = [
    GaussPoint {
        coords: C3D10_GAUSS_COORDS[0],
        weight: 0.041_666_666_666_667,
        n: &C3D10_N[0],
        dn: &C3D10_DN[0],
    },
    GaussPoint {
        coords: C3D10_GAUSS_COORDS[1],
        weight: 0.041_666_666_666_667,
        n: &C3D10_N[1],
        dn: &C3D10_DN[1],
    },
    GaussPoint {
        coords: C3D10_GAUSS_COORDS[2],
        weight: 0.041_666_666_666_667,
        n: &C3D10_N[2],
        dn: &C3D10_DN[2],
    },
    GaussPoint {
        coords: C3D10_GAUSS_COORDS[3],
        weight: 0.041_666_666_666_667,
        n: &C3D10_N[3],
        dn: &C3D10_DN[3],
    },
];

const C3D10_LOCAL_COORDS: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.0], // Node 1
    [1.0, 0.0, 0.0], // Node 2
    [0.0, 1.0, 0.0], // Node 3
    [0.0, 0.0, 1.0], // Node 4
    [0.5, 0.0, 0.0], // Node 5 (mid 1-2)
    [0.5, 0.5, 0.0], // Node 6 (mid 2-3)
    [0.0, 0.5, 0.0], // Node 7 (mid 3-1)
    [0.0, 0.0, 0.5], // Node 8 (mid 1-4)
    [0.5, 0.0, 0.5], // Node 9 (mid 2-4)
    [0.0, 0.5, 0.5], // Node 10 (mid 3-4)
];

/// Quadratic tetrahedron (C3D10).
#[derive(Debug, Clone)]
pub struct C3D10 {
    pub id: ElementId,
    pub nodes: [NodeId; 10],
}

impl FiniteElement for C3D10 {
    fn id(&self) -> ElementId {
        self.id
    }

    fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    fn num_nodes(&self) -> usize {
        10
    }

    fn vtk_cell_type(&self) -> u8 {
        24 // VTK_QUADRATIC_TETRA
    }

    fn integration_points(&self) -> &'static [GaussPoint] {
        &C3D10_GAUSS
    }

    fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (DVector<f64>, DMatrix<f64>) {
        let (n, dn) = c3d10_math(xi, eta, zeta);
        (
            DVector::from_column_slice(&n),
            DMatrix::from_column_slice(3, 10, &dn),
        )
    }

    fn node_local_coords(&self) -> &'static [[f64; 3]] {
        &C3D10_LOCAL_COORDS
    }
}

/// Mathematical definition of C3D10 shape functions and their derivatives.
///
/// Returns (N, dN) where N is a 10-element array and dN is a flattened 3x10 matrix
/// in COLUMN-MAJOR order.
const fn c3d10_math(xi: f64, et: f64, ze: f64) -> ([f64; 10], [f64; 30]) {
    let a = 1.0 - xi - et - ze;

    let n = [
        (2.0 * a - 1.0) * a,
        xi * (2.0 * xi - 1.0),
        et * (2.0 * et - 1.0),
        ze * (2.0 * ze - 1.0),
        4.0 * xi * a,
        4.0 * xi * et,
        4.0 * et * a,
        4.0 * ze * a,
        4.0 * xi * ze,
        4.0 * et * ze,
    ];

    let dn = [
        // Node 1
        1.0 - 4.0 * a, // dN/dxi
        1.0 - 4.0 * a, // dN/deta
        1.0 - 4.0 * a, // dN/dzeta
        // Node 2
        4.0 * xi - 1.0,
        0.0,
        0.0,
        // Node 3
        0.0,
        4.0 * et - 1.0,
        0.0,
        // Node 4
        0.0,
        0.0,
        4.0 * ze - 1.0,
        // Node 5
        4.0 * (a - xi),
        -4.0 * xi,
        -4.0 * xi,
        // Node 6
        4.0 * et,
        4.0 * xi,
        0.0,
        // Node 7
        -4.0 * et,
        4.0 * (a - et),
        -4.0 * et,
        // Node 8
        -4.0 * ze,
        -4.0 * ze,
        4.0 * (a - ze),
        // Node 9
        4.0 * ze,
        0.0,
        4.0 * xi,
        // Node 10
        0.0,
        4.0 * ze,
        4.0 * et,
    ];

    (n, dn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_func_sum() {
        let (n, _) = c3d10_math(0.25, 0.25, 0.25);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        let (n, _) = c3d10_math(0.1, 0.2, 0.3);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shape_func_nodes() {
        for (i, coord) in C3D10_LOCAL_COORDS.iter().enumerate() {
            let (n, _) = c3d10_math(coord[0], coord[1], coord[2]);
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
