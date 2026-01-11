use crate::solver::results::IncResult;
use std::error::Error;

/// Trait for writing output formats
pub trait ResultWriter {
    fn write_increment(&self, inc_result: &IncResult) -> Result<(), Box<dyn Error>>;
    /// Is called at the very end of the analysis. Can be used for cleanup, etc.
    fn finish(&self);
}
