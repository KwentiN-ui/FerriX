use derive_more::{Display, FromStr};

use crate::solver::io::{pvd::PvdWriter, writer::ResultWriter};

pub mod pvd;
pub mod vtk;
pub mod writer;

#[derive(Debug, Clone, Display, FromStr)]
pub enum OutputFormat {
    #[display("PVD")]
    Pvd,
}

impl OutputFormat {
    pub fn get_writer(&self) -> Box<dyn ResultWriter> {
        match self {
            OutputFormat::Pvd => Box::new(PvdWriter::new()),
        }
    }
}

