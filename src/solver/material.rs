//! This module contains data structures for material definitions. A Material is defined using the `*MATERIAL` card
//! globally. To use a material, it has to be referenced in a section definition.
#![allow(unused)]

use nalgebra::DMatrix;

use crate::solver::inp::InpFile;

#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub density: Option<f64>,
    /// Elastic Modulus, Poisson
    pub elastic: Option<(f64, f64)>,
}

impl Material {
    /// constructs linear elastic materialmatrix D (6x6 for 3D)
    /// Voigt-Notation: xx, yy, zz, xy, yz, zx
    pub fn build_elastic_d_matrix(&self) -> DMatrix<f64> {
        let (e, nu) = self.elastic.expect("No *ELASTIC Card definition found!");

        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c1 = 1.0 - nu;
        let c2 = nu;
        let c3 = (1.0 - 2.0 * nu) / 2.0;

        DMatrix::from_row_slice(
            6,
            6,
            &[
                c1, c2, c2, 0., 0., 0., c2, c1, c2, 0., 0., 0., c2, c2, c1, 0., 0., 0., 0., 0., 0.,
                c3, 0., 0., 0., 0., 0., 0., c3, 0., 0., 0., 0., 0., 0., c3,
            ],
        ) * factor
    }
}
