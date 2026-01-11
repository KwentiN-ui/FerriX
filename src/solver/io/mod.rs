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
    pub fn get_writer(&self, project: &Box<Project>) -> Box<dyn ResultWriter> {
        match self {
            OutputFormat::Vtk => Box::new(VtkWriter::new(project.clone())),
        }
    }
}
