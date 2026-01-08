// C3D20 Element
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
    let mut n = vec![0.0; 20];

    // Corner nodes (1-8)
    n[0] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta) * (-xi - eta - zeta - 2.0);
    n[1] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta) * (xi - eta - zeta - 2.0);
    n[2] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta) * (xi + eta - zeta - 2.0);
    n[3] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta) * (-xi + eta - zeta - 2.0);
    n[4] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta) * (-xi - eta + zeta - 2.0);
    n[5] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta) * (xi - eta + zeta - 2.0);
    n[6] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta) * (xi + eta + zeta - 2.0);
    n[7] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta) * (-xi + eta + zeta - 2.0);

    // Midside nodes (9-20)
    n[8] = 0.25 * (1.0 - xi * xi) * (1.0 - eta) * (1.0 - zeta);
    n[9] = 0.25 * (1.0 + xi) * (1.0 - eta * eta) * (1.0 - zeta);
    n[10] = 0.25 * (1.0 - xi * xi) * (1.0 + eta) * (1.0 - zeta);
    n[11] = 0.25 * (1.0 - xi) * (1.0 - eta * eta) * (1.0 - zeta);
    n[12] = 0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta * zeta);
    n[13] = 0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta * zeta);
    n[14] = 0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta * zeta);
    n[15] = 0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta * zeta);
    n[16] = 0.25 * (1.0 - xi * xi) * (1.0 - eta) * (1.0 + zeta);
    n[17] = 0.25 * (1.0 + xi) * (1.0 - eta * eta) * (1.0 + zeta);
    n[18] = 0.25 * (1.0 - xi * xi) * (1.0 + eta) * (1.0 + zeta);
    n[19] = 0.25 * (1.0 - xi) * (1.0 - eta * eta) * (1.0 + zeta);

    // Derivatives
    let mut dn = Array2::<f64>::zeros((3, 20));

    // dN/dxi
    dn[[0, 0]] = 0.125 * (1.0 - eta) * (1.0 - zeta) * (-(-xi - eta - zeta - 2.0) - (1.0 - xi));
    dn[[0, 1]] = 0.125 * (1.0 - eta) * (1.0 - zeta) * (1.0 * (xi - eta - zeta - 2.0) + (1.0 + xi));
    dn[[0, 2]] = 0.125 * (1.0 + eta) * (1.0 - zeta) * (1.0 * (xi + eta - zeta - 2.0) + (1.0 + xi));
    dn[[0, 3]] = 0.125 * (1.0 + eta) * (1.0 - zeta) * (-(-xi + eta - zeta - 2.0) - (1.0 - xi));
    dn[[0, 4]] = 0.125 * (1.0 - eta) * (1.0 + zeta) * (-(-xi - eta + zeta - 2.0) - (1.0 - xi));
    dn[[0, 5]] = 0.125 * (1.0 - eta) * (1.0 + zeta) * (1.0 * (xi - eta + zeta - 2.0) + (1.0 + xi));
    dn[[0, 6]] = 0.125 * (1.0 + eta) * (1.0 + zeta) * (1.0 * (xi + eta + zeta - 2.0) + (1.0 + xi));
    dn[[0, 7]] = 0.125 * (1.0 + eta) * (1.0 + zeta) * (-(-xi + eta + zeta - 2.0) - (1.0 - xi));
    dn[[0, 8]] = -0.5 * xi * (1.0 - eta) * (1.0 - zeta);
    dn[[0, 9]] = 0.25 * (1.0 - eta * eta) * (1.0 - zeta);
    dn[[0, 10]] = -0.5 * xi * (1.0 + eta) * (1.0 - zeta);
    dn[[0, 11]] = -0.25 * (1.0 - eta * eta) * (1.0 - zeta);
    dn[[0, 12]] = -0.25 * (1.0 - eta) * (1.0 - zeta * zeta);
    dn[[0, 13]] = 0.25 * (1.0 - eta) * (1.0 - zeta * zeta);
    dn[[0, 14]] = 0.25 * (1.0 + eta) * (1.0 - zeta * zeta);
    dn[[0, 15]] = -0.25 * (1.0 + eta) * (1.0 - zeta * zeta);
    dn[[0, 16]] = -0.5 * xi * (1.0 - eta) * (1.0 + zeta);
    dn[[0, 17]] = 0.25 * (1.0 - eta * eta) * (1.0 + zeta);
    dn[[0, 18]] = -0.5 * xi * (1.0 + eta) * (1.0 + zeta);
    dn[[0, 19]] = -0.25 * (1.0 - eta * eta) * (1.0 + zeta);

    // dN/deta
    dn[[1, 0]] = 0.125 * (1.0 - xi) * (1.0 - zeta) * (-(-xi - eta - zeta - 2.0) - (1.0 - eta));
    dn[[1, 1]] = 0.125 * (1.0 + xi) * (1.0 - zeta) * (-(xi - eta - zeta - 2.0) - (1.0 - eta));
    dn[[1, 2]] = 0.125 * (1.0 + xi) * (1.0 - zeta) * (1.0 * (xi + eta - zeta - 2.0) + (1.0 + eta));
    dn[[1, 3]] = 0.125 * (1.0 - xi) * (1.0 - zeta) * (1.0 * (-xi + eta - zeta - 2.0) + (1.0 + eta));
    dn[[1, 4]] = 0.125 * (1.0 - xi) * (1.0 + zeta) * (-(-xi - eta + zeta - 2.0) - (1.0 - eta));
    dn[[1, 5]] = 0.125 * (1.0 + xi) * (1.0 + zeta) * (-(xi - eta + zeta - 2.0) - (1.0 - eta));
    dn[[1, 6]] = 0.125 * (1.0 + xi) * (1.0 + zeta) * (1.0 * (xi + eta + zeta - 2.0) + (1.0 + eta));
    dn[[1, 7]] = 0.125 * (1.0 - xi) * (1.0 + zeta) * (1.0 * (-xi + eta + zeta - 2.0) + (1.0 + eta));
    dn[[1, 8]] = -0.25 * (1.0 - xi * xi) * (1.0 - zeta);
    dn[[1, 9]] = -0.5 * eta * (1.0 + xi) * (1.0 - zeta);
    dn[[1, 10]] = 0.25 * (1.0 - xi * xi) * (1.0 - zeta);
    dn[[1, 11]] = -0.5 * eta * (1.0 - xi) * (1.0 - zeta);
    dn[[1, 12]] = -0.25 * (1.0 - xi) * (1.0 - zeta * zeta);
    dn[[1, 13]] = -0.25 * (1.0 + xi) * (1.0 - zeta * zeta);
    dn[[1, 14]] = 0.25 * (1.0 + xi) * (1.0 - zeta * zeta);
    dn[[1, 15]] = 0.25 * (1.0 - xi) * (1.0 - zeta * zeta);
    dn[[1, 16]] = -0.25 * (1.0 - xi * xi) * (1.0 + zeta);
    dn[[1, 17]] = -0.5 * eta * (1.0 + xi) * (1.0 + zeta);
    dn[[1, 18]] = 0.25 * (1.0 - xi * xi) * (1.0 + zeta);
    dn[[1, 19]] = -0.5 * eta * (1.0 - xi) * (1.0 + zeta);

    // dN/dzeta
    dn[[2, 0]] = 0.125 * (1.0 - xi) * (1.0 - eta) * (-(-xi - eta - zeta - 2.0) - (1.0 - zeta));
    dn[[2, 1]] = 0.125 * (1.0 + xi) * (1.0 - eta) * (-(xi - eta - zeta - 2.0) - (1.0 - zeta));
    dn[[2, 2]] = 0.125 * (1.0 + xi) * (1.0 + eta) * (-(xi + eta - zeta - 2.0) - (1.0 - zeta));
    dn[[2, 3]] = 0.125 * (1.0 - xi) * (1.0 + eta) * (-(-xi + eta - zeta - 2.0) - (1.0 - zeta));
    dn[[2, 4]] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 * (-xi - eta + zeta - 2.0) + (1.0 + zeta));
    dn[[2, 5]] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 * (xi - eta + zeta - 2.0) + (1.0 + zeta));
    dn[[2, 6]] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 * (xi + eta + zeta - 2.0) + (1.0 + zeta));
    dn[[2, 7]] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 * (-xi + eta + zeta - 2.0) + (1.0 + zeta));
    dn[[2, 8]] = -0.25 * (1.0 - xi * xi) * (1.0 - eta);
    dn[[2, 9]] = -0.25 * (1.0 + xi) * (1.0 - eta * eta);
    dn[[2, 10]] = -0.25 * (1.0 - xi * xi) * (1.0 + eta);
    dn[[2, 11]] = -0.25 * (1.0 - xi) * (1.0 - eta * eta);
    dn[[2, 12]] = -0.5 * zeta * (1.0 - xi) * (1.0 - eta);
    dn[[2, 13]] = -0.5 * zeta * (1.0 + xi) * (1.0 - eta);
    dn[[2, 14]] = -0.5 * zeta * (1.0 + xi) * (1.0 + eta);
    dn[[2, 15]] = -0.5 * zeta * (1.0 - xi) * (1.0 + eta);
    dn[[2, 16]] = 0.25 * (1.0 - xi * xi) * (1.0 - eta);
    dn[[2, 17]] = 0.25 * (1.0 + xi) * (1.0 - eta * eta);
    dn[[2, 18]] = 0.25 * (1.0 - xi * xi) * (1.0 + eta);
    dn[[2, 19]] = 0.25 * (1.0 - xi) * (1.0 - eta * eta);

    (n, dn)
}
