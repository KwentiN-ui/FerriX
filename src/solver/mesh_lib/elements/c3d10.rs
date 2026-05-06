//! C3D10 element implementation (Quadratic Tetrahedron).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

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
        let (n, dn) = shape_func_c3d10(xi, eta, zeta);
        (
            DVector::from_column_slice(n.as_slice()),
            DMatrix::from_column_slice(3, 10, dn.as_slice()),
        )
    }

    fn node_local_coords(&self) -> &'static [[f64; 3]] {
        &C3D10_LOCAL_COORDS
    }
}

/// Shape functions for C3D10
#[must_use]
pub fn shape_func_c3d10(xi: f64, et: f64, ze: f64) -> (SVector<f64, 10>, SMatrix<f64, 3, 10>) {
    let a = 1.0 - xi - et - ze;

    let mut n = SVector::<f64, 10>::zeros();
    n[0] = (2.0 * a - 1.0) * a;
    n[1] = xi * (2.0 * xi - 1.0);
    n[2] = et * (2.0 * et - 1.0);
    n[3] = ze * (2.0 * ze - 1.0);
    n[4] = 4.0 * xi * a;
    n[5] = 4.0 * xi * et;
    n[6] = 4.0 * et * a;
    n[7] = 4.0 * ze * a;
    n[8] = 4.0 * xi * ze;
    n[9] = 4.0 * et * ze;

    let mut dn = SMatrix::<f64, 3, 10>::zeros();
    // xi-derivatives
    dn[(0, 0)] = 1.0 - 4.0 * a;
    dn[(0, 1)] = 4.0 * xi - 1.0;
    dn[(0, 2)] = 0.0;
    dn[(0, 3)] = 0.0;
    dn[(0, 4)] = 4.0 * (a - xi);
    dn[(0, 5)] = 4.0 * et;
    dn[(0, 6)] = -4.0 * et;
    dn[(0, 7)] = -4.0 * ze;
    dn[(0, 8)] = 4.0 * ze;
    dn[(0, 9)] = 0.0;

    // eta-derivatives
    dn[(1, 0)] = 1.0 - 4.0 * a;
    dn[(1, 1)] = 0.0;
    dn[(1, 2)] = 4.0 * et - 1.0;
    dn[(1, 3)] = 0.0;
    dn[(1, 4)] = -4.0 * xi;
    dn[(1, 5)] = 4.0 * xi;
    dn[(1, 6)] = 4.0 * (a - et);
    dn[(1, 7)] = -4.0 * ze;
    dn[(1, 8)] = 0.0;
    dn[(1, 9)] = 4.0 * ze;

    // zeta-derivatives
    dn[(2, 0)] = 1.0 - 4.0 * a;
    dn[(2, 1)] = 0.0;
    dn[(2, 2)] = 0.0;
    dn[(2, 3)] = 4.0 * ze - 1.0;
    dn[(2, 4)] = -4.0 * xi;
    dn[(2, 5)] = 0.0;
    dn[(2, 6)] = -4.0 * et;
    dn[(2, 7)] = 4.0 * (a - ze);
    dn[(2, 8)] = 4.0 * xi;
    dn[(2, 9)] = 4.0 * et;

    (n, dn)
}

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

/// Gaussian integration points for C3D10
const C3D10_GAUSS: [GaussPoint; 4] = [
    GaussPoint {
        coords: [
            0.138_196_601_125_011,
            0.138_196_601_125_011,
            0.138_196_601_125_011,
        ],
        weight: 0.041_666_666_666_667,
    },
    GaussPoint {
        coords: [
            0.585_410_196_624_968,
            0.138_196_601_125_011,
            0.138_196_601_125_011,
        ],
        weight: 0.041_666_666_666_667,
    },
    GaussPoint {
        coords: [
            0.138_196_601_125_011,
            0.585_410_196_624_968,
            0.138_196_601_125_011,
        ],
        weight: 0.041_666_666_666_667,
    },
    GaussPoint {
        coords: [
            0.138_196_601_125_011,
            0.138_196_601_125_011,
            0.585_410_196_624_968,
        ],
        weight: 0.041_666_666_666_667,
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_func_sum() {
        let (n, _) = shape_func_c3d10(0.25, 0.25, 0.25);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        let (n, _) = shape_func_c3d10(0.1, 0.2, 0.3);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shape_func_nodes() {
        let local_coords = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];

        for (i, coord) in local_coords.iter().enumerate() {
            let (n, _) = shape_func_c3d10(coord[0], coord[1], coord[2]);
            assert!(
                (n[i] - 1.0).abs() < 1e-12,
                "Node {} failed: n[{}] = {}",
                i + 1,
                i,
                n[i]
            );
            for j in 0..10 {
                if i != j {
                    assert!(
                        n[j].abs() < 1e-12,
                        "Node {} failed at n[{}]: n[{}] = {}",
                        i + 1,
                        j,
                        j,
                        n[j]
                    );
                }
            }
        }
    }
}
