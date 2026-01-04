use std::{error::Error, fs::read_to_string};

use ccx_rs::Step;

pub struct Project {
    steps: Vec<Box<dyn Step>>,
}

impl Project {
    pub fn from_filepath(path: &str) -> Result<Self, Box<dyn Error>> {
        let content = read_to_string(path)?;
        Self::from_str(&content)
    }

    pub fn from_str(string: &str) -> Result<Self, Box<dyn Error>> {
        todo!()
    }
}
