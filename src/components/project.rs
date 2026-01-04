use std::{error::Error, fs::read_to_string, path::Path};
use thiserror::Error;

use ccx_rs::Step;

use crate::components::mesh_lib::mesh::Mesh;

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
    pub fn from_str(string: &str) -> Result<Self, Box<dyn Error>> {
        let sections = parse_sections_from_str(string)?;
        println!("{sections:?}");
        todo!()
    }
}
fn parse_sections_from_str(string: &str) -> Result<Vec<InpSection>, Box<dyn Error>> {
    let mut sections: Vec<InpSection> = Vec::new();
    for (nr, line) in string.lines().enumerate() {
        // sanitize the line
        let line = sanitize_line(line);
        if line.starts_with("**") {
            // ignore comment
            continue;
        }
        if line == "*NODE" {
            sections.push(InpSection::Node(nr));
        } else if line == "*HEADING" {
            sections.push(InpSection::Heading(nr));
        }
    }
    check_if_valid(sections)
}

/// Checks the sections for obvious mistakes, like missing node definitions, etc...
fn check_if_valid(sections: Vec<InpSection>) -> Result<Vec<InpSection>, Box<dyn Error>> {
    // TODO: Should be expanded later
    // Check if the sections contain a Node Card
    if !sections
        .iter()
        .any(|sec| matches!(sec, InpSection::Node(_)))
    {
        return Err(Box::new(InpParsingError::MissingNodeDefinition));
    }
    Ok(sections)
}

fn sanitize_line(line: &str) -> String {
    line.trim()
        .chars()
        .map(|c| c.to_uppercase().to_string())
        .collect()
}

/// Different types of sections of the .inp file and their line-number
#[derive(Debug, PartialEq)]
pub enum InpSection {
    Node(usize),
    Heading(usize),
    Element(usize),
}

#[derive(Error, Debug)]
pub enum InpParsingError {
    #[error("Keyword {0} is unknown (Line {1})")]
    UnknownKeyword(String, u64),

    #[error("The input file does not contain a *NODE Card. Analysis is aborted.")]
    MissingNodeDefinition,
}
