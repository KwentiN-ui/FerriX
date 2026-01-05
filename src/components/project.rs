use std::{
    collections::HashMap,
    error::Error,
    fs::{self, read_to_string},
};
use thiserror::Error;

use crate::components::{
    mesh_lib::{
        elements::element::{Element, ElementType},
        mesh::Mesh,
    },
    step::step_trait::Step,
};

pub struct Project {
    pub mesh: Mesh,
    pub steps: Vec<Box<dyn Step>>,
}

impl Project {
    pub fn print_info(&self) {
        println!("--- Project Info ---");
        println!("Mesh:");
        println!("  Nodes: {}", self.mesh.nodes.len());
        println!("  Elements:");

        let mut elem_count: HashMap<ElementType, u32> = HashMap::new();
        for elem in &self.mesh.elements {
            let elem_type: ElementType = elem.into();
            *elem_count.entry(elem_type).or_insert(0) += 1;
        }
        for (elem, count) in elem_count {
            println!("  - {elem:?}: {count}");
        }
    }
    pub fn from_jobname(
        jobname: &str,
        preprocess_output: Option<&String>,
    ) -> Result<Self, Box<dyn Error>> {
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
        Self::from_str(&content, preprocess_output)
    }

    /// Attempts to parse an .inp file into a Project.
    pub fn from_str(
        raw_input: &str,
        preprocess_output: Option<&String>,
    ) -> Result<Self, Box<dyn Error>> {
        let input = preprocess_inp(raw_input);
        if let Some(out) = preprocess_output {
            fs::write(out, &input)?;
        }
        let sections = parse_sections_from_str(&input);

        let mesh = Mesh::from_sections(&input, &sections)?;
        Ok(Self {
            mesh,
            steps: Vec::new(),
        })
    }
}

fn parse_sections_from_str(string: &str) -> Vec<InpSection> {
    let mut sections: Vec<InpSection> = Vec::new();
    for (nr, line) in string.lines().enumerate() {
        if line == "*NODE" {
            sections.push(InpSection::Node(nr));
        } else if line.starts_with("*ELEMENT") {
            sections.push(InpSection::Element(nr));
        }
    }
    sections
}

/// This preprocess includes:
/// - removing leading and trailing whitespaces
/// - removing comments
/// - making all text uppercase
/// - merging lines that belong together (`,` at the end of line)
/// - removes empty lines
fn preprocess_inp(input_file: &str) -> String {
    input_file
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(|line| line.chars().map(|c| c.to_uppercase().to_string()).collect())
        // merge lines that end with `,`
        .map(|line: String| {
            if line.ends_with(',') {
                line + " "
            } else {
                line + "\n"
            }
        })
        // remove comments
        .filter(|line| !line.starts_with("**"))
        .collect()
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
        let inp = "**comment\n word  \n \t*keyword\n123.4\n4, 5, 6,\n7, 8, 9";
        assert_eq!(
            preprocess_inp(inp),
            "WORD\n*KEYWORD\n123.4\n4, 5, 6, 7, 8, 9\n"
        );
    }
}
