use crate::solver::{results::IncResult, time::SolverTime};
use std::error::Error;

/// Trait for writing output formats
pub trait ResultWriter {
    /// Is called at the beginning of analysis. Can be used to setup directories etc.
    ///
    /// # Errors
    /// Returns an error if the initialization fails (e.g. directory creation).
    fn init(&self) -> Result<(), Box<dyn Error>>;

    /// Writes the results of an increment.
    ///
    /// # Errors
    /// Returns an error if the writing fails.
    fn write_increment(
        &self,
        inc_result: &IncResult,
        timer: &SolverTime,
    ) -> Result<(), Box<dyn Error>>;

    /// Is called at the very end of the analysis. Can be used for cleanup, etc.
    ///
    /// # Errors
    /// Returns an error if the finish operation fails.
    fn finish(&self) -> Result<(), Box<dyn Error>>;
}
