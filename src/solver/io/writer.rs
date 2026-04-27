//! Trait definitions for result exporters.

use crate::solver::{results::IncResult, time::SolverTime};
use std::error::Error;

/// A generic interface for writing simulation results to a file or stream.
pub trait ResultWriter {
    /// Initializes the writer at the start of the analysis.
    ///
    /// This can be used to create output directories, write file headers,
    /// or initialize data structures.
    ///
    /// # Errors
    /// Returns an error if the initialization fails (e.g., directory creation).
    fn init(&self) -> Result<(), Box<dyn Error>>;

    /// Writes the results of a single completed simulation increment.
    ///
    /// # Errors
    /// Returns an error if writing to the output medium fails.
    fn write_increment(
        &self,
        inc_result: &IncResult,
        timer: &SolverTime,
    ) -> Result<(), Box<dyn Error>>;

    /// Finalizes the output at the end of the analysis.
    ///
    /// This can be used to close files, write footers, or perform cleanup.
    ///
    /// # Errors
    /// Returns an error if the finalization fails.
    fn finish(&self) -> Result<(), Box<dyn Error>>;
}
