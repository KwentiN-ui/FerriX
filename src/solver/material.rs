//! Material property definitions.
//!
//! This module contains structures for defining physical properties like density
//! and elasticity, and provides methods for generating material law matrices.

use nalgebra::DMatrix;

/// Defines the physical and mechanical properties of a material.
#[derive(Debug, Clone)]
pub struct Material {
    /// Unique name of the material.
    pub name: String,
    /// Mass density of the material.
    pub density: Option<f64>,
    /// Elastic properties: (Young's modulus E, Poisson's ratio nu).
    pub elastic: Option<(f64, f64)>,
}

impl Material {
    /// Builds the elastic constitutive matrix (D-matrix) for the material (6x6 for 3D).
    /// Uses Voigt notation: [xx, yy, zz, xy, yz, zx].
    ///
    /// # Panics
    /// Panics if the material does not have an elastic definition (`elastic` is `None`).
    #[must_use]
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
                c1, c2, c2, 0.0, 0.0, 0.0, c2, c1, c2, 0.0, 0.0, 0.0, c2, c2, c1, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, c3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, c3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                c3,
            ],
        ) * factor
    }
}
