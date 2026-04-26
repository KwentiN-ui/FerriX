use crate::solver::error::{FerrixError, Result};
use crate::solver::solvers::Solver;
use faer::Side;
use faer::prelude::*;
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};

pub struct DirectSolver;

impl Solver for DirectSolver {
    fn solve(
        &self,
        k_global: &sprs::CsMat<f64>,
        f_global: &[f64],
        _preconditioner: Option<&dyn crate::solver::preconditioner::Preconditioner>,
        _tol: f64,
        _max_iter: usize,
    ) -> Result<Vec<f64>> {
        let (rows, cols) = k_global.shape();

        let indptr_storage = k_global.indptr();
        let indptr = indptr_storage.as_slice().ok_or_else(|| {
            FerrixError::NumericalError("Failed to get indptr from sprs matrix".into())
        })?;

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
            FerrixError::NumericalError(format!(
                "Cholesky factorization failed: {e}. Check for rigid body modes!"
            ))
        })?;

        let x_mat = llt.solve(&b);

        Ok(x_mat.col(0).iter().copied().collect())
    }
}
