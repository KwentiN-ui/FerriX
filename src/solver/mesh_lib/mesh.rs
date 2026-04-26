use crate::solver::ids::{ElementId, NodeId};
use std::collections::HashMap;

use crate::solver::{
    mesh_lib::{
        elements::element::{Element, ElementType},
        node::Node,
    },
};

/// Contains all Node and Element Data
#[derive(Debug, Clone, Default)]
pub struct Mesh {
    pub nodes: HashMap<NodeId, Node>,
    pub elements: HashMap<ElementId, Element>,
    pub node_sets: HashMap<String, Vec<NodeId>>,
    pub element_sets: HashMap<String, Vec<ElementId>>,

    /// Map: Node-ID (from INP) -> Matrix-Index (0..N)
    pub node_id_to_index: HashMap<NodeId, usize>,
    /// Map: Matrix-Index (0..N) -> Node-ID (from INP)
    pub index_to_node_id: Vec<NodeId>,
}

impl Mesh {
    pub fn build_node_mappings(&mut self) {
        let mut sorted_ids: Vec<NodeId> = self.nodes.keys().copied().collect();
        // Cannot derive Ord for NodeId, so we sort by the inner value
        sorted_ids.sort_unstable_by_key(|a| a.0);

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
    #[must_use] 
    pub fn get_index_for_node_id(&self, id: NodeId) -> Option<usize> {
        self.node_id_to_index.get(&id).copied()
    }

    /// Counts all elements by their respective type.
    #[must_use] 
    pub fn count_by_type(&self) -> HashMap<ElementType, u32> {
        let mut elem_count: HashMap<ElementType, u32> = HashMap::new();
        for elem in self.elements.values() {
            let elem_type: ElementType = elem.into();
            *elem_count.entry(elem_type).or_insert(0) += 1;
        }
        elem_count
    }
}
