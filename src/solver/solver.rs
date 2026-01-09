use crate::solver::preconditioner::Preconditioner;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::ArrayView1;
use rayon::prelude::*;
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
    /// * `preconditioner` - An optional preconditioner to accelerate convergence.
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
        preconditioner: Option<&dyn Preconditioner>,
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String>;
}

/// An iterative solver using the Preconditioned Conjugate Gradient (PCG) method.
///
/// This is a specific implementation of the `Solver` trait that uses the PCG method,
/// which is well-suited for large, sparse, symmetric, and positive-definite systems
/// commonly found in FEA.
pub struct IterativeSolver;

impl Solver for IterativeSolver {
    fn solve(
        &self,
        k: &CsMat<f64>,
        b: &[f64],
        preconditioner: Option<&dyn Preconditioner>,
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String> {
        println!("Conjugate gradient solver is running...");

        let spinner = ProgressBar::new_spinner();
        spinner.enable_steady_tick(Duration::from_millis(100));
        spinner.set_style(ProgressStyle::default_spinner());

        let b_len = b.len();
        let mut x = vec![0.0; b_len]; // Initial guess x0 = 0

        let mut r = b.to_vec();
        let mut z = if let Some(p) = preconditioner {
            p.apply(&r)
        } else {
            r.clone()
        };
        let mut p = z.clone();
        
        let mut rs_old: f64 = r.par_iter().zip(&z).map(|(&r_val, &z_val)| r_val * z_val).sum();

        for i in 0..max_iter {
            spinner.set_message(format!("Iteration: {i}, Residual: {rs_old:.6}"));
            if rs_old.sqrt() < tol {
                spinner.finish();
                return Ok(x);
            }

            let p_view = ArrayView1::from(&p);
            let ap = (k * &p_view).to_vec();

            let alpha = rs_old / p.par_iter().zip(&ap).map(|(&p_val, &ap_val)| p_val * ap_val).sum::<f64>();

            x.par_iter_mut().zip(&p).for_each(|(x_val, &p_val)| {
                *x_val += alpha * p_val;
            });
            r.par_iter_mut().zip(&ap).for_each(|(r_val, &ap_val)| {
                *r_val -= alpha * ap_val;
            });

            if let Some(precond) = preconditioner {
                z = precond.apply(&r);
            } else {
                z = r.clone();
            }

            let rs_new: f64 = r.par_iter().zip(&z).map(|(&r_val, &z_val)| r_val * z_val).sum();
            
            let beta = rs_new / rs_old;

            p.par_iter_mut().zip(&z).for_each(|(p_val, &z_val)| {
                *p_val = z_val + beta * *p_val;
            });

            rs_old = rs_new;
        }
        spinner.finish();
        Err(format!(
            "CG solver did not converge after {max_iter} iterations"
        ))
    }
}
