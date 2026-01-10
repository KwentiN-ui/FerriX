use std::time::Duration;

use crate::solver::solvers::Solver;
use faer::Side;
use faer::prelude::*;
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

pub struct DirectSolver;

impl Solver for DirectSolver {
    fn solve(
        &self,
        k_global: &sprs::CsMat<f64>,
        f_global: &[f64],
        _preconditioner: Option<&dyn crate::solver::preconditioner::Preconditioner>,
        _tol: f64,
        _max_iter: usize,
    ) -> Result<Vec<f64>, String> {
        let spinner = ProgressBar::new_spinner();
        spinner.set_message("Solving using Direct Solver (Cholesky-Decomposition)");
        spinner.enable_steady_tick(Duration::from_millis(100));
        spinner.set_style(ProgressStyle::default_spinner());

        let (rows, cols) = k_global.shape();

        let indptr_storage = k_global.indptr();
        let indptr = indptr_storage
            .as_slice()
            .ok_or("Failed to get indptr from sprs matrix")?;

        let indices = k_global.indices();
        let data = k_global.data();

        let k_faer = SparseColMatRef::<'_, usize, f64>::new(
            SymbolicSparseColMatRef::new_checked(rows, cols, indptr, None, indices),
            data,
        );

        // prepare right side (F)
        let b = faer::Mat::<f64>::from_fn(f_global.len(), 1, |i, _| f_global[i]);

        // Cholesky-Decomposition
        let llt = k_faer.sp_cholesky(Side::Lower).map_err(|e| {
            format!("Cholesky factorization failed: {e}. Check for rigid body modes!")
        })?;

        let x_mat = llt.solve(&b);
        spinner.finish();

        Ok(x_mat.col(0).iter().copied().collect())
    }
}
