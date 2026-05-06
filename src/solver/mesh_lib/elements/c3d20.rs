//! C3D20 element implementation (Quadratic Hexahedron).

use crate::solver::ids::{ElementId, NodeId};
use crate::solver::mesh_lib::elements::element::{FiniteElement, GaussPoint};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

const C3D20_GAUSS: [GaussPoint; 27] = compute_gauss();

/// Compiletime creation of Gauss-Points
const fn compute_gauss() -> [GaussPoint; 27] {
    let pts = [-0.774_596_669_241_483, 0.0, 0.774_596_669_241_483];
    let wts = [
        0.555_555_555_555_555_6,
        0.888_888_888_888_888_8,
        0.555_555_555_555_555_6,
    ];

    // Initialisiere ein Array mit Platzhalter-Werten
    let mut gps = [GaussPoint {
        coords: [0.0; 3],
        weight: 0.0,
    }; 27];
    let mut index = 0;

    // while-Schleifen sind in const-Kontexten problemlos möglich
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut k = 0;
            while k < 3 {
                gps[index] = GaussPoint {
                    coords: [pts[i], pts[j], pts[k]],
                    weight: wts[i] * wts[j] * wts[k],
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
        let (n, dn) = shape_func_c3d20(xi, et, ze);
        (
            DVector::from_column_slice(n.as_slice()),
            DMatrix::from_column_slice(3, 20, dn.as_slice()),
        )
    }

    fn node_local_coords(&self) -> &'static [[f64; 3]] {
        &C3D20_LOCAL_COORDS
    }
}

/// Legacy function for `compute_stiffness_sdv` glue code
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn shape_func_c3d20(xi: f64, et: f64, ze: f64) -> (SVector<f64, 20>, SMatrix<f64, 3, 20>) {
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

    let omgmopg = (omg - opg) / 4.0;
    let omhmoph = (omh - oph) / 4.0;
    let omrmopr = (omr - opr) / 4.0;

    let mut n = SVector::<f64, 20>::zeros();

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

    let mut dn = SMatrix::<f64, 3, 20>::zeros();

    // xi-derivatives
    dn[(0, 0)] = omh * omr * (tpgphpr - omg) / 8.0;
    dn[(0, 1)] = (opg - tmgphpr) * omh * omr / 8.0;
    dn[(0, 2)] = (opg - tmgmhpr) * oph * omr / 8.0;
    dn[(0, 3)] = oph * omr * (tpgmhpr - omg) / 8.0;
    dn[(0, 4)] = omh * opr * (tpgphmr - omg) / 8.0;
    dn[(0, 5)] = (opg - tmgphmr) * omh * opr / 8.0;
    dn[(0, 6)] = (opg - tmgmhmr) * oph * opr / 8.0;
    dn[(0, 7)] = oph * opr * (tpgmhmr - omg) / 8.0;
    dn[(0, 8)] = omgmopg * omh * omr;
    dn[(0, 9)] = omhoph * omr;
    dn[(0, 10)] = omgmopg * oph * omr;
    dn[(0, 11)] = -omhoph * omr;
    dn[(0, 12)] = omgmopg * omh * opr;
    dn[(0, 13)] = omhoph * opr;
    dn[(0, 14)] = omgmopg * oph * opr;
    dn[(0, 15)] = -omhoph * opr;
    dn[(0, 16)] = -omropr * omh;
    dn[(0, 17)] = omropr * omh;
    dn[(0, 18)] = omropr * oph;
    dn[(0, 19)] = -omropr * oph;

    // eta-derivatives
    dn[(1, 0)] = omg * omr * (tpgphpr - omh) / 8.0;
    dn[(1, 1)] = opg * omr * (tmgphpr - omh) / 8.0;
    dn[(1, 2)] = opg * (oph - tmgmhpr) * omr / 8.0;
    dn[(1, 3)] = omg * (oph - tpgmhpr) * omr / 8.0;
    dn[(1, 4)] = omg * opr * (tpgphmr - omh) / 8.0;
    dn[(1, 5)] = opg * opr * (tmgphmr - omh) / 8.0;
    dn[(1, 6)] = opg * (oph - tmgmhmr) * opr / 8.0;
    dn[(1, 7)] = omg * (oph - tpgmhmr) * opr / 8.0;
    dn[(1, 9)] = omhmoph * opg * omr;
    dn[(1, 8)] = -omgopg * omr;
    dn[(1, 10)] = omgopg * omr;
    dn[(1, 11)] = omhmoph * omg * omr;
    dn[(1, 12)] = -omgopg * opr;
    dn[(1, 13)] = omhmoph * opg * opr;
    dn[(1, 14)] = omgopg * opr;
    dn[(1, 15)] = omhmoph * omg * opr;
    dn[(1, 16)] = -omropr * omg;
    dn[(1, 17)] = -omropr * opg;
    dn[(1, 18)] = omropr * opg;
    dn[(1, 19)] = omropr * omg;

    // zeta-derivatives
    dn[(2, 0)] = omg * omh * (tpgphpr - omr) / 8.0;
    dn[(2, 1)] = opg * omh * (tmgphpr - omr) / 8.0;
    dn[(2, 2)] = opg * oph * (tmgmhpr - omr) / 8.0;
    dn[(2, 3)] = omg * oph * (tpgmhpr - omr) / 8.0;
    dn[(2, 4)] = omg * omh * (opr - tpgphmr) / 8.0;
    dn[(2, 5)] = opg * omh * (opr - tmgphmr) / 8.0;
    dn[(2, 6)] = opg * oph * (opr - tmgmhmr) / 8.0;
    dn[(2, 7)] = omg * oph * (opr - tpgmhmr) / 8.0;
    dn[(2, 8)] = -omgopg * omh;
    dn[(2, 9)] = -omhoph * opg;
    dn[(2, 10)] = -omgopg * oph;
    dn[(2, 11)] = -omhoph * omg;
    dn[(2, 12)] = omgopg * omh;
    dn[(2, 13)] = omhoph * opg;
    dn[(2, 14)] = omgopg * oph;
    dn[(2, 15)] = omhoph * omg;
    dn[(2, 16)] = omrmopr * omg * omh;
    dn[(2, 17)] = omrmopr * opg * omh;
    dn[(2, 18)] = omrmopr * opg * oph;
    dn[(2, 19)] = omrmopr * omg * oph;

    (n, dn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_func_sum() {
        let (n, _) = shape_func_c3d20(0.0, 0.0, 0.0);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        let (n, _) = shape_func_c3d20(0.5, -0.3, 0.8);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shape_func_corners() {
        // At corner 1 (-1, -1, -1) -> node 1 index 0
        let (n, _) = shape_func_c3d20(-1.0, -1.0, -1.0);
        assert!((n[0] - 1.0).abs() < 1e-12);
        for i in 1..20 {
            assert!(n[i].abs() < 1e-12);
        }

        // At mid-side 9 (0, -1, -1) -> index 8
        let (n, _) = shape_func_c3d20(0.0, -1.0, -1.0);
        assert!((n[8] - 1.0).abs() < 1e-12);
    }
}
