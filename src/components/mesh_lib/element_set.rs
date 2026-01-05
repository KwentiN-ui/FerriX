use std::error::Error;

use crate::components::mesh_lib::elements::element::Element;

#[derive(Debug)]
pub struct ElementSet {
    elem_type: String,
    elements: Vec<Element>,
}

impl ElementSet {
    pub fn from_string(string: &str, line_nr: usize) -> Result<Self, Box<dyn Error>> {
        let elem_type: String = string
            .lines()
            .nth(line_nr)
            .expect("This should not be able to happen. The line_nr is outside the file-range")
            .split(',')
            .map(str::trim)
            .nth(1)
            .ok_or("Invalid Element definition on line {line_nr}")?
            .split('=')
            .next_back()
            .ok_or("Invalid Element definition on line {line_nr}")?
            .to_string();

        let elements = Vec::new();
        Ok(Self {
            elem_type,
            elements,
        })
    }
}
