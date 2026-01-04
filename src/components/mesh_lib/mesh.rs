use std::collections::HashMap;

use crate::components::mesh_lib::{element_set::ElementSet, node::Node};

/// Contains all Node and Element Data
pub struct Mesh {
    nodes: HashMap<u64, Node>,
    element_sets: Vec<ElementSet>,
}
