//! A trait for preconditioners used in iterative solvers.

use rayon::prelude::*;
use sprs::CsMat;

/// A trait for preconditioners.
///
/// A preconditioner `M` is a matrix that approximates the inverse of the global
/// stiffness matrix `K`. It is used to transform the system of linear equations
/// into a form that is easier to solve, thus accelerating convergence.
pub trait Preconditioner {
    /// Applies the preconditioner to a vector.
    ///
    /// # Arguments
    ///
    /// * `r` - The residual vector to be preconditioned.
    ///
    /// # Returns
    ///
    /// The preconditioned residual vector.
    fn apply(&self, r: &[f64]) -> Vec<f64>;
}

/// A simple diagonal (or Jacobi) preconditioner.
///
/// This preconditioner uses the inverse of the diagonal of the stiffness matrix `K`.
/// It is computationally inexpensive and can be effective for some problems.
pub struct DiagonalPreconditioner {
    inv_diag: Vec<f64>,
}

impl DiagonalPreconditioner {
    /// Creates a new `DiagonalPreconditioner` from the global stiffness matrix.
    #[must_use] 
    pub fn new(k_global: &CsMat<f64>) -> Self {
        let n = k_global.rows();
        let mut inv_diag = vec![1.0; n]; // Default 1.0 für leere Zeilen

        for (i, &val) in k_global.diag().iter() {
            if val.abs() > 1e-9 {
                inv_diag[i] = 1.0 / val;
            }
        }
        Self { inv_diag }
    }
}

impl Preconditioner for DiagonalPreconditioner {
    fn apply(&self, r: &[f64]) -> Vec<f64> {
        r.par_iter()
            .zip(&self.inv_diag)
            .map(|(&r_val, &inv_d_val)| r_val * inv_d_val)
            .collect()
    }
}
