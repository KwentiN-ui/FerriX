// Placeholder for the C3D20 Element

#![allow(unused)]

use ndarray::Array2;

use crate::solver::mesh_lib::elements::element::GaussPoint;

pub fn c3d20_gauss() -> Vec<GaussPoint> {
    // Gauss-Legendre
    let val = 1.0 / (3.0_f64).sqrt();
    let w = 1.0;
    let mut points = Vec::with_capacity(8);
    for k in [-val, val] {
        for j in [-val, val] {
            for i in [-val, val] {
                points.push(GaussPoint {
                    coords: [i, j, k],
                    weight: w * w * w,
                });
            }
        }
    }
    points
}

pub fn shape_func_c3d20(xi: f64, eta: f64, zeta: f64) -> (Vec<f64>, Array2<f64>) {
    todo!()
}
