use std::{error::Error, fs::read_to_string, path::Path};
use thiserror::Error;

use crate::components::{mesh_lib::mesh::Mesh, step::step_trait::Step};

pub struct Project {
    mesh: Mesh,
    steps: Vec<Box<dyn Step>>,
}

impl Project {
    pub fn from_jobname(jobname: &str) -> Result<Self, Box<dyn Error>> {
        // check if the jobname already has a file extension
        let path = if std::path::Path::new(jobname)
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("inp"))
        {
            jobname.to_string()
        } else {
            jobname.to_string() + ".inp"
        };
        let content = read_to_string(path)?;
        Self::from_str(&content)
    }

    /// Attempts to parse an .inp file into a Project.
    pub fn from_str(raw_input: &str) -> Result<Self, Box<dyn Error>> {
        let input = preprocess_inp(raw_input);
        let sections = parse_sections_from_str(&input);

        let mesh = Mesh::from_sections(&input, &sections);
        todo!()
    }
}

fn parse_sections_from_str(string: &str) -> Vec<InpSection> {
    let mut sections: Vec<InpSection> = Vec::new();
    for (nr, line) in string.lines().enumerate() {
        // sanitize the line
        if line.starts_with("**") {
            // ignore comment
            continue;
        }
        if line == "*NODE" {
            sections.push(InpSection::Node(nr));
        } else if line.starts_with("*ELEMENT") {
            sections.push(InpSection::Element(nr));
        }
    }
    sections
}

/// Removes whitespaces and comments, and makes all text uppercase
fn preprocess_inp(input_file: &str) -> String {
    input_file
        .lines()
        .map(str::trim)
        .map(make_allcaps)
        .map(|line| line + "\n")
        // remove comments
        .filter(|line| !line.starts_with("**"))
        .collect()
}

fn make_allcaps(line: &str) -> String {
    line.chars().map(|c| c.to_uppercase().to_string()).collect()
}

/// Different types of sections of the .inp file and their line-number
#[derive(Debug, PartialEq)]
pub enum InpSection {
    Node(usize),
    Element(usize),
}

#[derive(Error, Debug)]
pub enum InpParsingError {
    #[error("Keyword {0} is unknown (Line {1})")]
    UnknownKeyword(String, u64),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_preprocess_inp() {
        let inp = "**comment\n word  \n \t*keyword\n123.4";
        assert_eq!(preprocess_inp(inp), "WORD\n*KEYWORD\n123.4\n");
    }
}
