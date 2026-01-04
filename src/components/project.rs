use std::{error::Error, fs::read_to_string};

use ccx_rs::Step;

use crate::components::mesh_lib::mesh::Mesh;

pub struct Project {
    mesh: Mesh,
    steps: Vec<Box<dyn Step>>,
}

impl Project {
    pub fn from_filepath(path: &str) -> Result<Self, Box<dyn Error>> {
        let content = read_to_string(path)?;
        Self::from_str(&content)
    }

    /// Attempts to parse an .inp file into a Project.
    pub fn from_str(string: &str) -> Result<Self, Box<dyn Error>> {
        todo!()
    }
}
