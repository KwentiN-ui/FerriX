use ndarray::{Array2, array};

use crate::solver::mesh_lib::elements::element::GaussPoint;

/// linear tetraeder
pub fn shape_func_c3d4(xi: f64, eta: f64, zeta: f64) -> (Vec<f64>, Array2<f64>) {
    // Shape functions (CalculiX Convention: Node 1 is Origin)
    // N1 = 1 - xi - eta - zeta
    // N2 = xi
    // N3 = eta
    // N4 = zeta
    let n = vec![
        1.0 - xi - eta - zeta, // Node 1
        xi,                    // Node 2
        eta,                   // Node 3
        zeta,                  // Node 4
    ];

    // Derivatives
    // Row 0: d/dxi
    // Row 1: d/deta
    // Row 2: d/dzeta
    //
    // cols are nodes 1, 2, 3, 4
    let dn = array![
        [-1.0, 1.0, 0.0, 0.0], // d/dxi:  N1=-1, N2=1
        [-1.0, 0.0, 1.0, 0.0], // d/deta: N1=-1, N3=1
        [-1.0, 0.0, 0.0, 1.0]  // d/dzeta: N1=-1, N4=1
    ];

    (n, dn)
}

pub fn c3d4_gauss() -> Vec<GaussPoint> {
    vec![
        GaussPoint {
            coords: [0.25, 0.25, 0.25],
            weight: 1.0 / 6.0,
        }, // Volumen Tetraeder = 1/6
    ]
}
