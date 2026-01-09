use super::{Solver, Preconditioner};
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::ArrayView1;
use rayon::prelude::*;
use sprs::CsMat;
use std::time::Duration;

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

        let mut residual = b.to_vec();
        let mut preconditioned_residual = if let Some(p) = preconditioner {
            p.apply(&residual)
        } else {
            residual.clone()
        };
        let mut search_direction = preconditioned_residual.clone();
        
        let mut rs_old: f64 = residual.par_iter().zip(&preconditioned_residual).map(|(&r_val, &z_val)| r_val * z_val).sum();

        for i in 0..max_iter {
            spinner.set_message(format!("Iteration: {i}, Residual: {rs_old:.6}"));
            if rs_old.sqrt() < tol {
                spinner.finish();
                return Ok(x);
            }

            let search_direction_view = ArrayView1::from(&search_direction);
            let ap = (k * &search_direction_view).to_vec();

            let alpha = rs_old / search_direction.par_iter().zip(&ap).map(|(&p_val, &ap_val)| p_val * ap_val).sum::<f64>();

            x.par_iter_mut().zip(&search_direction).for_each(|(x_val, &p_val)| {
                *x_val += alpha * p_val;
            });
            residual.par_iter_mut().zip(&ap).for_each(|(r_val, &ap_val)| {
                *r_val -= alpha * ap_val;
            });

            if let Some(precond) = preconditioner {
                preconditioned_residual = precond.apply(&residual);
            } else {
                preconditioned_residual.clone_from(&residual);
            }

            let rs_new: f64 = residual.par_iter().zip(&preconditioned_residual).map(|(&r_val, &z_val)| r_val * z_val).sum();
            
            let beta = rs_new / rs_old;

            search_direction.par_iter_mut().zip(&preconditioned_residual).for_each(|(p_val, &z_val)| {
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
