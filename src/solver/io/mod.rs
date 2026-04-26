use std::sync::Arc;

use derive_more::{Display, FromStr};

use crate::solver::{
    io::{vtk::VtkWriter, writer::ResultWriter},
    project::Project,
};

pub mod vtk;
pub mod writer;

#[derive(Debug, Clone, Display, FromStr)]
pub enum OutputFormat {
    #[display("VTK")]
    Vtk,
}

impl OutputFormat {
    #[must_use] 
    pub fn get_writer(&self, project: Arc<Project>) -> Box<dyn ResultWriter> {
        match self {
            OutputFormat::Vtk => Box::new(VtkWriter::new(project)),
        }
    }
}
