//! C3D4 element implementation (Linear Tetrahedron).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

const C3D4_GAUSS: [GaussPoint; 1] = [GaussPoint {
    coords: [0.25, 0.25, 0.25],
    weight: 1.0 / 6.0,
}];

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
        let n = SVector::<f64, 4>::new(1.0 - xi - eta - zeta, xi, eta, zeta);

        let dn = SMatrix::<f64, 3, 4>::from_row_slice(&[
            -1.0, 1.0, 0.0, 0.0, // Row 0: dN/dxi
            -1.0, 0.0, 1.0, 0.0, // Row 1: dN/deta
            -1.0, 0.0, 0.0, 1.0, // Row 2: dN/dzeta
        ]);

        (
            DVector::from_column_slice(n.as_slice()),
            DMatrix::from_column_slice(3, 4, dn.as_slice()),
        )
    }

    fn node_local_coords(&self) -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
}

/// Legacy function for `compute_stiffness_sdv` glue code
#[must_use]
pub fn shape_func_c3d4(xi: f64, eta: f64, zeta: f64) -> (SVector<f64, 4>, SMatrix<f64, 3, 4>) {
    let n = SVector::<f64, 4>::new(1.0 - xi - eta - zeta, xi, eta, zeta);

    let dn = SMatrix::<f64, 3, 4>::from_row_slice(&[
        -1.0, 1.0, 0.0, 0.0, // Row 0: dN/dxi
        -1.0, 0.0, 1.0, 0.0, // Row 1: dN/deta
        -1.0, 0.0, 0.0, 1.0, // Row 2: dN/dzeta
    ]);

    (n, dn)
}

/// Legacy function for `compute_stiffness_sdv` glue code
#[must_use]
pub fn c3d4_gauss() -> Vec<GaussPoint> {
    vec![GaussPoint {
        coords: [0.25, 0.25, 0.25],
        weight: 1.0 / 6.0,
    }]
}
