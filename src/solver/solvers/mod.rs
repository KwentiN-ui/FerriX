//! Linear system solvers.
//!
//! This module defines the `Solver` trait and provides implementations for
//! direct (Cholesky) and iterative (Preconditioned Conjugate Gradient) methods.

pub mod direct;
pub mod iterative;

use crate::solver::error::Result;
use crate::solver::preconditioner::Preconditioner;
use sprs::CsMat;

/// A trait for linear system solvers.
///
/// This trait defines a generic interface for solving a system of linear equations
/// of the form `K * u = F`, where `K` is the global stiffness matrix, `u` is the
/// displacement vector to be solved for, and `F` is the global force vector.
pub trait Solver {
    /// Solves the system of linear equations.
    ///
    /// # Arguments
    ///
    /// * `k_global` - The global stiffness matrix `K` as a compressed sparse row (CSR) matrix.
    /// * `f_global` - The global force vector `F` as a slice.
    /// * `preconditioner` - An optional preconditioner to accelerate convergence.
    /// * `tol` - The tolerance for convergence.
    /// * `max_iter` - The maximum number of iterations allowed.
    ///
    /// # Returns
    ///
    /// A `Result` containing the displacement vector `u` if the solution converges,
    /// or a `FerrixError` otherwise.
    ///
    /// # Errors
    /// Returns an error if the solver fails to converge or encounters a numerical issue.
    fn solve(
        &self,
        k_global: &CsMat<f64>,
        f_global: &[f64],
        preconditioner: Option<&dyn Preconditioner>,
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>>;
}

/// Strategies for solving the global system of equations.
#[derive(Debug, Clone)]
pub enum SolverType {
    /// A direct solver using sparse Cholesky factorization.
    Direct,
    /// An iterative solver using the Preconditioned Conjugate Gradient (PCG) method.
    Iterative,
    /// The default solver strategy for the current analysis type.
    Default,
}
