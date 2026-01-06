use std::{collections::HashMap, error::Error};

use crate::solver::{
    inp::InpFile,
    mesh_lib::{
        elements::element::{Element, ElementType},
        node::Node,
    },
    project::InpSection,
};

/// Contains all Node and Element Data
#[derive(Debug, Clone)]
pub struct Mesh {
    pub nodes: HashMap<usize, Node>,
    pub elements: HashMap<usize, Element>,
    pub node_sets: HashMap<String, Vec<usize>>,

    /// Map: Node-ID (from INP) -> Matrix-Index (0..N)
    pub node_id_to_index: HashMap<usize, usize>,
    /// Map: Matrix-Index (0..N) -> Node-ID (from INP)
    pub index_to_node_id: Vec<usize>,
}

impl Mesh {
    pub fn from_sections(
        input_file: &InpFile,
        sections: &Vec<InpSection>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut elements = Vec::new();
        let mut nodes: Option<Vec<Node>> = None;
        let mut node_sets: HashMap<String, Vec<usize>> = HashMap::new();

        for sec in sections {
            match sec {
                InpSection::Node(nr) => {
                    nodes = Some(
                        input_file
                            .0
                            .lines()
                            .skip(*nr + 1)
                            .map_while(Node::parse_line)
                            .collect(),
                    );
                }
                InpSection::Element(nr) => {
                    let elem_type = Element::parse_type_str_from_line(
                        input_file
                            .0
                            .lines()
                            .nth(*nr)
                            .expect("Line number outside file"),
                    )?;
                    elements.extend(
                        input_file
                            .0
                            .lines()
                            .skip(nr + 1)
                            .take_while(|line| {
                                line.chars().nth(0).expect("No empty lines").is_numeric()
                            })
                            .map(|line| Element::parse_line(&elem_type, line)),
                    );
                }
                // Parse Node Sets
                // Assumes `InpSection::Nset(line_index, name, is_generate)`
                InpSection::Nset(nr, name, is_generate) => {
                    let mut ids = Vec::new();

                    // Read lines until next keyword (*)
                    let data_lines = input_file
                        .0
                        .lines()
                        .skip(nr + 1)
                        .take_while(|line| !line.trim().starts_with('*'));

                    if *is_generate {
                        // Format: Start, End, Increment
                        for line in data_lines {
                            let nums: Vec<usize> = line
                                .split(',')
                                .map(str::trim)
                                .filter(|s| !s.is_empty())
                                .filter_map(|s| s.parse().ok())
                                .collect();

                            if nums.len() >= 2 {
                                // At least Start, End
                                let start = nums[0];
                                let end = nums[1];
                                let step = *nums.get(2).unwrap_or(&1); // Default step 1

                                for id in (start..=end).step_by(step) {
                                    ids.push(id);
                                }
                            }
                        }
                    } else {
                        // Explicit List: 1, 2, 3, 4...
                        for line in data_lines {
                            let line_ids = line
                                .split(',')
                                .map(str::trim)
                                .filter(|s| !s.is_empty())
                                .filter_map(|s| s.parse::<usize>().ok());
                            ids.extend(line_ids);
                        }
                    }
                    node_sets.insert(name.clone(), ids);
                }
            }
        }

        let mut node_hash: HashMap<usize, Node> = HashMap::new();
        for node in nodes.ok_or("Input file missing *NODE card.")? {
            node_hash.insert(node.id, node);
        }

        let mut elem_hash: HashMap<usize, Element> = HashMap::new();
        for elem in elements {
            elem_hash.insert(elem.get_id(), elem);
        }

        // 3. Construct Mesh with populated sets
        let mut mesh = Self {
            nodes: node_hash,
            elements: elem_hash,
            node_sets, // Pass the filled map
            node_id_to_index: HashMap::new(),
            index_to_node_id: Vec::new(),
        };

        mesh.build_node_mappings();

        Ok(mesh)
    }

    pub fn build_node_mappings(&mut self) {
        let mut sorted_ids: Vec<usize> = self.nodes.keys().copied().collect();
        sorted_ids.sort_unstable();

        // Reverse Mapping (Index -> ID)
        self.index_to_node_id.clone_from(&sorted_ids);

        // Forward Mapping (ID -> Index)
        self.node_id_to_index = sorted_ids
            .into_iter()
            .enumerate()
            .map(|(idx, id)| (id, idx))
            .collect();
    }

    /// Fetches the matrix-index for a given Node-ID (from INP file)
    pub fn get_index_for_node_id(&self, id: usize) -> Option<usize> {
        self.node_id_to_index.get(&id).copied()
    }

    /// Counts all elements by their respective type.
    pub fn count_by_type(&self) -> HashMap<ElementType, u32> {
        let mut elem_count: HashMap<ElementType, u32> = HashMap::new();
        for elem in self.elements.values() {
            let elem_type: ElementType = elem.into();
            *elem_count.entry(elem_type).or_insert(0) += 1;
        }
        elem_count
    }
}
