use indicatif::{ProgressBar, ProgressStyle};
use ndarray::ArrayView1;
use sprs::CsMat;
use std::time::Duration;

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
    /// * `tol` - The tolerance for convergence.
    /// * `max_iter` - The maximum number of iterations allowed.
    ///
    /// # Returns
    ///
    /// A `Result` containing the displacement vector `u` if the solution converges,
    /// or a `String` with an error message otherwise.
    fn solve(
        &self,
        k_global: &CsMat<f64>,
        f_global: &[f64],
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String>;
}

/// An iterative solver using the Conjugate Gradient (CG) method.
///
/// This is a specific implementation of the `Solver` trait that uses the CG method,
/// which is well-suited for large, sparse, symmetric, and positive-definite systems
/// commonly found in FEA.
pub struct IterativeSolver;

impl Solver for IterativeSolver {
    fn solve(
        &self,
        k: &CsMat<f64>,
        b: &[f64],
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String> {
        println!("Conjugate gradient solver is running...");

        let spinner = ProgressBar::new_spinner();
        spinner.enable_steady_tick(Duration::from_millis(100));
        spinner.set_style(ProgressStyle::default_spinner());

        let b_len = b.len();
        let mut x = vec![0.0; b_len]; // Initial guess x0 = 0

        // r = b - K * x0  (since x0=0 -> r=b)
        let mut r = b.to_vec();
        let mut p = r.clone();

        // r_old = r^T * r
        let mut rs_old: f64 = r.iter().map(|v| v * v).sum();

        for _ in 0..max_iter {
            spinner.set_message(format!("Residual: {rs_old:.6}"));
            if rs_old.sqrt() < tol {
                spinner.finish();
                return Ok(x);
            }

            // Ap = K * p
            let p_view = ArrayView1::from(&p);
            let ap = (k * &p_view).to_vec();

            // alpha = rs_old / (p^T * Ap)
            let p_dot_ap: f64 = p.iter().zip(&ap).map(|(pi, api)| pi * api).sum();
            if p_dot_ap.abs() < 1e-15 {
                return Err("CG breakdown: denominator zero (matrix singular?)".to_string());
            }
            let alpha = rs_old / p_dot_ap;

            // x = x + alpha * p
            // r = r - alpha * Ap
            for j in 0..b_len {
                x[j] += alpha * p[j];
                r[j] -= alpha * ap[j];
            }

            let rs_new: f64 = r.iter().map(|v| v * v).sum();

            // beta = rs_new / rs_old
            let beta = rs_new / rs_old;

            // p = r + beta * p
            for j in 0..b_len {
                p[j] = r[j] + beta * p[j];
            }

            rs_old = rs_new;
        }
        spinner.finish();
        Err(format!(
            "CG solver did not converge after {max_iter} iterations"
        ))
    }
}
