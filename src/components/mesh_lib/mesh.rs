use std::{collections::HashMap, error::Error};

use crate::components::{
    mesh_lib::{elements::element::Element, node::Node},
    project::{InpParsingError, InpSection},
};

/// Contains all Node and Element Data
#[derive(Debug)]
pub struct Mesh {
    pub nodes: HashMap<usize, Node>,
    pub elements: Vec<Element>,
}

impl Mesh {
    #[allow(clippy::match_wildcard_for_single_variants)]
    pub fn from_sections(
        input_file: &str,
        sections: &Vec<InpSection>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut elements = Vec::new();
        let mut nodes: Option<Vec<Node>> = None;
        for sec in sections {
            match sec {
                InpSection::Node(nr) => {
                    nodes = Some(
                        input_file
                            .lines()
                            .skip(*nr + 1)
                            .map_while(Node::parse_line)
                            .collect(),
                    );
                }
                InpSection::Element(nr) => {
                    let elem_type = Element::parse_type_str_from_line(
                        input_file
                            .lines()
                            .nth(*nr)
                            .expect("The line number is outside the file, aborting"),
                    )?;
                    elements.extend(
                        input_file
                            .lines()
                            .skip(nr + 1)
                            .take_while(|line| {
                                line.chars()
                                    .nth(0)
                                    .expect("There are no empty lines after preprocessing")
                                    .is_numeric()
                            })
                            .map(|line| Element::parse_line(&elem_type, line)),
                    );
                }

                _ => {}
            }
        }
        let mut node_hash: HashMap<usize, Node> = HashMap::new();
        for node in
            nodes.ok_or("The input file does not contain a *NODE card. Analysis is aborted.")?
        {
            node_hash.insert(node.id, node);
        }

        Ok(Self {
            nodes: node_hash,
            elements,
        })
    }
}
