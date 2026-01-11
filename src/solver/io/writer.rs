use crate::solver::project::Project;
use crate::solver::results::StepResult;
use std::error::Error;

/// Trait for writing output formats
pub trait ResultWriter {
    /// Called once before the main solver loop.
    fn init(&mut self, project: &Project) -> Result<(), Box<dyn Error>>;

    /// Called for each increment.
    fn write(&mut self, result: &StepResult) -> Result<(), Box<dyn Error>>;

    /// Called once after the main solver loop.
    fn finish(&mut self) -> Result<(), Box<dyn Error>>;
}
