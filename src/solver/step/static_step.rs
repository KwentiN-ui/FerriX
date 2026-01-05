use std::{error::Error, sync::Arc};

use crate::solver::{inp::InpFile, mesh_lib::mesh::Mesh, step::steps::Step};

#[derive(Debug, Clone)]
pub struct StaticStep {
    input: Arc<InpFile>,
    mesh: Arc<Mesh>,
}

impl StaticStep {
    pub fn new(input: Arc<InpFile>, mesh: Arc<Mesh>) -> Self {
        Self { input, mesh }
    }
    pub fn compute(&mut self) -> Result<(), Box<dyn Error>> {
        todo!()
    }
}
