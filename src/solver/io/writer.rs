use crate::solver::{results::IncResult, time::SolverTime};
use std::error::Error;

/// Trait for writing output formats
pub trait ResultWriter {
    /// Is called at the beginning of analysis. Can be used to setup directories etc.
    fn init(&self) -> Result<(), Box<dyn Error>>;
    fn write_increment(
        &self,
        inc_result: &IncResult,
        timer: &SolverTime,
    ) -> Result<(), Box<dyn Error>>;
    /// Is called at the very end of the analysis. Can be used for cleanup, etc.
    fn finish(&self) -> Result<(), Box<dyn Error>>;
}
