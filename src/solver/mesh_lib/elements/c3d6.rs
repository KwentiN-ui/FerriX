//! C3D6 element implementation (Linear Wedge).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

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
        let (n, dn) = shape_func_c3d6(xi, eta, zeta);
        (
            DVector::from_column_slice(n.as_slice()),
            DMatrix::from_column_slice(3, 6, dn.as_slice()),
        )
    }

    fn node_local_coords(&self) -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, -1.0], // Node 1
            [1.0, 0.0, -1.0], // Node 2
            [0.0, 1.0, -1.0], // Node 3
            [0.0, 0.0, 1.0],  // Node 4
            [1.0, 0.0, 1.0],  // Node 5
            [0.0, 1.0, 1.0],  // Node 6
        ]
    }
}

/// Shape functions for C3D6
#[must_use]
pub fn shape_func_c3d6(xi: f64, et: f64, ze: f64) -> (SVector<f64, 6>, SMatrix<f64, 3, 6>) {
    let a = 1.0 - xi - et;

    let mut n = SVector::<f64, 6>::zeros();
    n[0] = 0.5 * a * (1.0 - ze);
    n[1] = 0.5 * xi * (1.0 - ze);
    n[2] = 0.5 * et * (1.0 - ze);
    n[3] = 0.5 * a * (1.0 + ze);
    n[4] = 0.5 * xi * (1.0 + ze);
    n[5] = 0.5 * et * (1.0 + ze);

    let mut dn = SMatrix::<f64, 3, 6>::zeros();
    // xi-derivatives
    dn[(0, 0)] = -0.5 * (1.0 - ze);
    dn[(0, 1)] = 0.5 * (1.0 - ze);
    dn[(0, 2)] = 0.0;
    dn[(0, 3)] = -0.5 * (1.0 + ze);
    dn[(0, 4)] = 0.5 * (1.0 + ze);
    dn[(0, 5)] = 0.0;

    // eta-derivatives
    dn[(1, 0)] = -0.5 * (1.0 - ze);
    dn[(1, 1)] = 0.0;
    dn[(1, 2)] = 0.5 * (1.0 - ze);
    dn[(1, 3)] = -0.5 * (1.0 + ze);
    dn[(1, 4)] = 0.0;
    dn[(1, 5)] = 0.5 * (1.0 + ze);

    // zeta-derivatives
    dn[(2, 0)] = -0.5 * a;
    dn[(2, 1)] = -0.5 * xi;
    dn[(2, 2)] = -0.5 * et;
    dn[(2, 3)] = 0.5 * a;
    dn[(2, 4)] = 0.5 * xi;
    dn[(2, 5)] = 0.5 * et;

    (n, dn)
}

/// Gaussian integration points for C3D6
const C3D6_GAUSS: [GaussPoint; 2] = [
    GaussPoint {
        coords: [
            0.333_333_333_333_333,
            0.333_333_333_333_333,
            -0.577_350_269_189_626,
        ],
        weight: 0.5,
    },
    GaussPoint {
        coords: [
            0.333_333_333_333_333,
            0.333_333_333_333_333,
            0.577_350_269_189_626,
        ],
        weight: 0.5,
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_func_sum() {
        let (n, _) = shape_func_c3d6(0.25, 0.25, 0.0);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        let (n, _) = shape_func_c3d6(0.1, 0.2, 0.5);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shape_func_nodes() {
        let local_coords = vec![
            [0.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ];

        for (i, coord) in local_coords.iter().enumerate() {
            let (n, _) = shape_func_c3d6(coord[0], coord[1], coord[2]);
            assert!(
                (n[i] - 1.0).abs() < 1e-12,
                "Node {} failed: n[{}] = {}",
                i + 1,
                i,
                n[i]
            );
            for j in 0..6 {
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
