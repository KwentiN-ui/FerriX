use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::results::StepResult;
use std::error::Error;
use std::path::Path;

/// Trait for writing output formats
pub trait ResultWriter {
    fn write(
        &self,
        path: &Path,
        mesh: &Mesh,
        step_result: &[StepResult],
        nodal_output: &[String],
        element_output: &[String],
    ) -> Result<(), Box<dyn Error>>;
}
