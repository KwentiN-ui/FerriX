use std::{collections::HashMap, error::Error};

use crate::components::{
    mesh_lib::{element_set::ElementSet, node::Node},
    project::{InpParsingError, InpSection},
};

/// Contains all Node and Element Data
#[derive(Debug)]
pub struct Mesh {
    nodes: HashMap<usize, Node>,
    element_sets: Vec<ElementSet>,
}

impl Mesh {
    #[allow(clippy::match_wildcard_for_single_variants)]
    pub fn from_sections(
        input_file: &str,
        sections: &Vec<InpSection>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut element_sets = Vec::new();
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
                    element_sets.push(ElementSet::from_string(input_file, *nr)?);
                }

                _ => {}
            }
        }
        // TODO
        let mut node_hash: HashMap<usize, Node> = HashMap::new();
        for node in
            nodes.ok_or("The input file does not contain a *NODE card. Analysis is aborted.")?
        {
            node_hash.insert(node.id, node);
        }

        Ok(Self {
            nodes: node_hash,
            element_sets,
        })
    }
}
