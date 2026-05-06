//! C3D20 element implementation (Quadratic Hexahedron).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector};

/// Static storage for precomputed shape function values and derivatives for all 27 integration points.
static C3D20_N: [[f64; 20]; 27] = compute_all_n();
static C3D20_DN: [[f64; 60]; 27] = compute_all_dn();

/// Compile-time precomputation of shape function values (N) for all 27 integration points.
const fn compute_all_n() -> [[f64; 20]; 27] {
    let pts = [-0.774_596_669_241_483, 0.0, 0.774_596_669_241_483];
    let mut all_n = [[0.0; 20]; 27];
    let mut index = 0;
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut k = 0;
            while k < 3 {
                all_n[index] = c3d20_math(pts[i], pts[j], pts[k]).0;
                index += 1;
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    all_n
}

/// Compile-time precomputation of shape function derivatives (dN) for all 27 integration points.
const fn compute_all_dn() -> [[f64; 60]; 27] {
    let pts = [-0.774_596_669_241_483, 0.0, 0.774_596_669_241_483];
    let mut all_dn = [[0.0; 60]; 27];
    let mut index = 0;
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut k = 0;
            while k < 3 {
                all_dn[index] = c3d20_math(pts[i], pts[j], pts[k]).1;
                index += 1;
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    all_dn
}

/// Compile-time creation of Gauss-Points.
const fn compute_gauss() -> [GaussPoint; 27] {
    let pts = [-0.774_596_669_241_483, 0.0, 0.774_596_669_241_483];
    let wts = [
        0.555_555_555_555_555_6,
        0.888_888_888_888_888_8,
        0.555_555_555_555_555_6,
    ];

    let mut gps = [GaussPoint {
        coords: [0.0; 3],
        weight: 0.0,
        n: &C3D20_N[0],
        dn: &C3D20_DN[0],
    }; 27];
    let mut index = 0;

    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut k = 0;
            while k < 3 {
                gps[index] = GaussPoint {
                    coords: [pts[i], pts[j], pts[k]],
                    weight: wts[i] * wts[j] * wts[k],
                    n: &C3D20_N[index],
                    dn: &C3D20_DN[index],
                };
                index += 1;
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }

    gps
}

const C3D20_GAUSS: [GaussPoint; 27] = compute_gauss();

const C3D20_LOCAL_COORDS: [[f64; 3]; 20] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
    [0.0, -1.0, -1.0],
    [1.0, 0.0, -1.0],
    [0.0, 1.0, -1.0],
    [-1.0, 0.0, -1.0],
    [0.0, -1.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [-1.0, 0.0, 1.0],
    [-1.0, -1.0, 0.0],
    [1.0, -1.0, 0.0],
    [1.0, 1.0, 0.0],
    [-1.0, 1.0, 0.0],
];

/// Quadratic hexahedron (C3D20).
#[derive(Debug, Clone)]
pub struct C3D20 {
    pub id: ElementId,
    pub nodes: [NodeId; 20],
}

impl FiniteElement for C3D20 {
    fn id(&self) -> ElementId {
        self.id
    }

    fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    fn num_nodes(&self) -> usize {
        20
    }

    fn vtk_cell_type(&self) -> u8 {
        25 // VTK_QUADRATIC_HEXAHEDRON
    }

    fn integration_points(&self) -> &'static [GaussPoint] {
        &C3D20_GAUSS
    }

    fn shape_functions(&self, xi: f64, et: f64, ze: f64) -> (DVector<f64>, DMatrix<f64>) {
        let (n, dn) = c3d20_math(xi, et, ze);
        (
            DVector::from_column_slice(&n),
            DMatrix::from_column_slice(3, 20, &dn),
        )
    }

    fn node_local_coords(&self) -> &'static [[f64; 3]] {
        &C3D20_LOCAL_COORDS
    }
}

/// Mathematical definition of C3D20 shape functions and their derivatives.
///
/// Returns (N, dN) where N is a 20-element array and dN is a flattened 3x20 matrix
/// in COLUMN-MAJOR order.
#[must_use]
#[allow(clippy::too_many_lines)]
pub const fn c3d20_math(xi: f64, et: f64, ze: f64) -> ([f64; 20], [f64; 60]) {
    let omg = 1.0 - xi;
    let omh = 1.0 - et;
    let omr = 1.0 - ze;
    let opg = 1.0 + xi;
    let oph = 1.0 + et;
    let opr = 1.0 + ze;

    let tpgphpr = opg + oph + ze;
    let tmgphpr = omg + oph + ze;
    let tmgmhpr = omg + omh + ze;
    let tpgmhpr = opg + omh + ze;

    let tpgphmr = opg + oph - ze;
    let tmgphmr = omg + oph - ze;
    let tmgmhmr = omg + omh - ze;
    let tpgmhmr = opg + omh - ze;

    let omgopg = omg * opg / 4.0;
    let omhoph = omh * oph / 4.0;
    let omropr = omr * opr / 4.0;

    let mut n = [0.0; 20];

    // Corner nodes
    n[0] = -omg * omh * omr * tpgphpr / 8.0;
    n[1] = -opg * omh * omr * tmgphpr / 8.0;
    n[2] = -opg * oph * omr * tmgmhpr / 8.0;
    n[3] = -omg * oph * omr * tpgmhpr / 8.0;
    n[4] = -omg * omh * opr * tpgphmr / 8.0;
    n[5] = -opg * omh * opr * tmgphmr / 8.0;
    n[6] = -opg * oph * opr * tmgmhmr / 8.0;
    n[7] = -omg * oph * opr * tpgmhmr / 8.0;

    // Mid-side nodes
    n[8] = omgopg * omh * omr;
    n[9] = omhoph * opg * omr;
    n[10] = omgopg * oph * omr;
    n[11] = omhoph * omg * omr;

    n[12] = omgopg * omh * opr;
    n[13] = omhoph * opg * opr;
    n[14] = omgopg * oph * opr;
    n[15] = omhoph * omg * opr;

    n[16] = omropr * omg * omh;
    n[17] = omropr * opg * omh;
    n[18] = omropr * opg * oph;
    n[19] = omropr * omg * oph;

    let mut dn = [0.0; 60];

    // Node 1
    dn[0] = omh * omr * (tpgphpr - omg) / 8.0;
    dn[1] = omg * omr * (tpgphpr - omh) / 8.0;
    dn[2] = omg * omh * (tpgphpr - omr) / 8.0;

    // Node 2
    dn[3] = (opg - tmgphpr) * omh * omr / 8.0;
    dn[4] = opg * omr * (tmgphpr - omh) / 8.0;
    dn[5] = opg * omh * (tmgphpr - omr) / 8.0;

    // Node 3
    dn[6] = (opg - tmgmhpr) * oph * omr / 8.0;
    dn[7] = opg * (oph - tmgmhpr) * omr / 8.0;
    dn[8] = opg * oph * (tmgmhpr - omr) / 8.0;

    // Node 4
    dn[9] = oph * omr * (tpgmhpr - omg) / 8.0;
    dn[10] = omg * (oph - tpgmhpr) * omr / 8.0;
    dn[11] = omg * oph * (tpgmhpr - omr) / 8.0;

    // Node 5
    dn[12] = omh * opr * (tpgphmr - omg) / 8.0;
    dn[13] = omg * opr * (tpgphmr - omh) / 8.0;
    dn[14] = omg * omh * (opr - tpgphmr) / 8.0;

    // Node 6
    dn[15] = (opg - tmgphmr) * omh * opr / 8.0;
    dn[16] = opg * opr * (tmgphmr - omh) / 8.0;
    dn[17] = opg * omh * (opr - tmgphmr) / 8.0;

    // Node 7
    dn[18] = (opg - tmgmhmr) * oph * opr / 8.0;
    dn[19] = opg * (oph - tmgmhmr) * opr / 8.0;
    dn[20] = opg * oph * (opr - tmgmhmr) / 8.0;

    // Node 8
    dn[21] = oph * opr * (tpgmhmr - omg) / 8.0;
    dn[22] = omg * (oph - tpgmhmr) * opr / 8.0;
    dn[23] = omg * oph * (opr - tpgmhmr) / 8.0;

    // Mid-side xi-derivatives
    let omgmopg = (omg - opg) / 4.0;
    dn[24] = omgmopg * omh * omr; // Node 9 dN/dxi
    dn[25] = -omgopg * omr; // Node 9 dN/deta
    dn[26] = -omgopg * omh; // Node 9 dN/dzeta

    dn[27] = omhoph * omr; // Node 10 dN/dxi
    let omhmoph = (omh - oph) / 4.0;
    dn[28] = omhmoph * opg * omr; // Node 10 dN/deta
    dn[29] = -omhoph * opg; // Node 10 dN/dzeta

    dn[30] = omgmopg * oph * omr; // Node 11 dN/dxi
    dn[31] = omgopg * omr; // Node 11 dN/deta
    dn[32] = -omgopg * oph; // Node 11 dN/dzeta

    dn[33] = -omhoph * omr; // Node 12 dN/dxi
    dn[34] = omhmoph * omg * omr; // Node 12 dN/deta
    dn[35] = -omhoph * omg; // Node 12 dN/dzeta

    dn[36] = omgmopg * omh * opr; // Node 13 dN/dxi
    dn[37] = -omgopg * opr; // Node 13 dN/deta
    dn[38] = omgopg * omh; // Node 13 dN/dzeta

    dn[39] = omhoph * opr; // Node 14 dN/dxi
    dn[40] = omhmoph * opg * opr; // Node 14 dN/deta
    dn[41] = omhoph * opg; // Node 14 dN/dzeta

    dn[42] = omgmopg * oph * opr; // Node 15 dN/dxi
    dn[43] = omgopg * opr; // Node 15 dN/deta
    dn[44] = omgopg * oph; // Node 15 dN/dzeta

    dn[45] = -omhoph * opr; // Node 16 dN/dxi
    dn[46] = omhmoph * omg * opr; // Node 16 dN/deta
    dn[47] = omhoph * omg; // Node 16 dN/dzeta

    let omrmopr = (omr - opr) / 4.0;
    dn[48] = -omropr * omh; // Node 17 dN/dxi
    dn[49] = -omropr * omg; // Node 17 dN/deta
    dn[50] = omrmopr * omg * omh; // Node 17 dN/dzeta

    dn[51] = omropr * omh; // Node 18 dN/dxi
    dn[52] = -omropr * opg; // Node 18 dN/deta
    dn[53] = omrmopr * opg * omh; // Node 18 dN/dzeta

    dn[54] = omropr * oph; // Node 19 dN/dxi
    dn[55] = omropr * opg; // Node 19 dN/deta
    dn[56] = omrmopr * opg * oph; // Node 19 dN/dzeta

    dn[57] = -omropr * oph; // Node 20 dN/dxi
    dn[58] = omropr * omg; // Node 20 dN/deta
    dn[59] = omrmopr * omg * oph; // Node 20 dN/dzeta

    (n, dn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_func_sum() {
        let (n, _) = c3d20_math(0.0, 0.0, 0.0);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        let (n, _) = c3d20_math(0.5, -0.3, 0.8);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shape_func_corners() {
        // At corner 1 (-1, -1, -1) -> node 1 index 0
        let (n, _) = c3d20_math(-1.0, -1.0, -1.0);
        assert!((n[0] - 1.0).abs() < 1e-12);
        for ni in n.iter().skip(1) {
            assert!(ni.abs() < 1e-12);
        }

        // At mid-side 9 (0, -1, -1) -> index 8
        let (n, _) = c3d20_math(0.0, -1.0, -1.0);
        assert!((n[8] - 1.0).abs() < 1e-12);
    }
}
