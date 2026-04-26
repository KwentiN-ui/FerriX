use super::{Preconditioner, Solver};
use crate::solver::error::{FerrixError, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use sprs::CsMat;
use std::time::Duration;

pub struct IterativeSolver;

impl Solver for IterativeSolver {
    fn solve(
        &self,
        k: &CsMat<f64>,
        b: &[f64],
        preconditioner: Option<&dyn Preconditioner>,
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>> {
        let spinner = ProgressBar::new_spinner();
        spinner.enable_steady_tick(Duration::from_millis(100));
        spinner.set_style(ProgressStyle::default_spinner());

        let b_len = b.len();
        let mut x = vec![0.0; b_len]; // Initial guess u0 = 0
        let mut residual = b.to_vec(); // r = b - K*x0

        // 1. Initial Preconditioning: z0 = M^-1 * r0
        let mut preconditioned_residual = if let Some(p) = preconditioner {
            p.apply(&residual)
        } else {
            residual.clone()
        };

        // Initial Search Direction: p0 = z0
        let mut search_direction = preconditioned_residual.clone();

        // Scalar product for convergence and beta calculation: rs_old = r * z
        let mut rs_old: f64 = residual
            .par_iter()
            .zip(&preconditioned_residual)
            .map(|(&r, &z)| r * z)
            .sum();

        for i in 0..max_iter {
            spinner.set_message(format!("Iter: {i}, Res: {:.3e}", rs_old.sqrt()));

            // Convergence check: ||r|| < tol
            if rs_old.sqrt() < tol {
                spinner.finish_with_message("Converged");
                return Ok(x);
            }

            // 2. Matrix-Vector Product: ap = K * p
            let mut ap = vec![0.0; b_len];
            multiply_sparse_dense(k, &search_direction, &mut ap);

            // 3. Step size alpha: alpha = (r * z) / (p * K * p)
            let p_kp: f64 = search_direction
                .par_iter()
                .zip(&ap)
                .map(|(&p, &akp)| p * akp)
                .sum();

            let alpha = rs_old / p_kp;

            // 4. Update Solution and Residual
            // x = x + alpha * p
            // r = r - alpha * K * p
            x.par_iter_mut()
                .zip(&search_direction)
                .for_each(|(xi, &pi)| *xi += alpha * pi);
            residual
                .par_iter_mut()
                .zip(&ap)
                .for_each(|(ri, &api)| *ri -= alpha * api);

            // 5. Apply Preconditioner to new residual: z_new = M^-1 * r_new
            if let Some(precond) = preconditioner {
                preconditioned_residual = precond.apply(&residual);
            } else {
                preconditioned_residual.clone_from(&residual);
            }

            // 6. Beta for next search direction
            let rs_new: f64 = residual
                .par_iter()
                .zip(&preconditioned_residual)
                .map(|(&r, &z)| r * z)
                .sum();

            let beta = rs_new / rs_old;

            // 7. Update Search Direction: p_new = z_new + beta * p_old
            search_direction
                .par_iter_mut()
                .zip(&preconditioned_residual)
                .for_each(|(pi, &zi)| *pi = zi + beta * *pi);

            rs_old = rs_new;
        }

        spinner.finish();
        Err(FerrixError::ConvergenceError(format!(
            "PCG did not converge after {max_iter} iterations"
        )))
    }
}

/// Performs y = A * x where A is a sparse matrix in CSR format and x/y are dense vectors.
fn multiply_sparse_dense(a: &CsMat<f64>, x: &[f64], y: &mut [f64]) {
    y.par_iter_mut().enumerate().for_each(|(row_idx, y_val)| {
        let mut sum = 0.0;
        // Direct acces to next row of CSR matrix
        if let Some(row) = a.outer_view(row_idx) {
            for (col_idx, &val) in row.iter() {
                sum += val * x[col_idx];
            }
        }
        *y_val = sum;
    });
}
