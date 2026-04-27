//! Input/Output operations and result exporters.
//!
//! This module handles the export of simulation results to various formats
//! (like VTK) and provides a generic interface for result writers.

use derive_more::{Display, FromStr};
use std::sync::Arc;

use crate::solver::{
    io::{vtk::VtkWriter, writer::ResultWriter},
    project::Project,
};

pub mod vtk;
pub mod writer;

/// Supported formats for exporting simulation results.
#[derive(Debug, Clone, Display, FromStr)]
pub enum OutputFormat {
    /// Visualization Toolkit (VTK) format, compatible with `ParaView`.
    #[display("VTK")]
    Vtk,
}

impl OutputFormat {
    /// Returns the appropriate result writer for the format.
    #[must_use]
    pub fn get_writer(&self, project: Arc<Project>) -> Box<dyn ResultWriter> {
        match self {
            OutputFormat::Vtk => Box::new(VtkWriter::new(project)),
        }
    }
}
