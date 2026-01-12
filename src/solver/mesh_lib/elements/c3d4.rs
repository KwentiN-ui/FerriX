use crate::solver::mesh_lib::elements::element::GaussPoint;
use nalgebra::{SMatrix, SVector};

/// Linear tetrahedron (C3D4)
/// Returns shape functions N and their derivatives dN/d(xi, eta, zeta)
///
/// The shape functions are:
/// $N_1 = 1 - \xi - \eta - \zeta$
/// $N_2 = \xi$
/// $N_3 = \eta$
/// $N_4 = \zeta$
pub fn shape_func_c3d4(xi: f64, eta: f64, zeta: f64) -> (SVector<f64, 4>, SMatrix<f64, 3, 4>) {
    // SVector and SMatrix are stack-allocated
    let n = SVector::<f64, 4>::new(1.0 - xi - eta - zeta, xi, eta, zeta);

    let dn = SMatrix::<f64, 3, 4>::from_row_slice(&[
        -1.0, 1.0, 0.0, 0.0, // Row 0: dN/dxi
        -1.0, 0.0, 1.0, 0.0, // Row 1: dN/deta
        -1.0, 0.0, 0.0, 1.0, // Row 2: dN/dzeta
    ]);

    (n, dn)
}

/// Returns the integration points for a C3D4 element.
/// For a linear tetrahedron, a single point at the centroid is sufficient.
pub fn c3d4_gauss() -> Vec<GaussPoint> {
    vec![GaussPoint {
        coords: [0.25, 0.25, 0.25],
        weight: 1.0 / 6.0, // Reference volume of a unit tetrahedron
    }]
}
