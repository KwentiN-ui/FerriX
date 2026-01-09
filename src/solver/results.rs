use crate::solver::ids::NodeId;
use std::collections::HashMap;

/// Types of fields we can store
#[allow(unused)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldType {
    Displacement, // U
    Stress,       // S
    Strain,       // E
}

/// Holds values for one specific field (e.g. Displacements for all nodes)
#[derive(Debug, Clone)]
pub struct NodalResult {
    #[allow(dead_code)]
    pub name: String, // e.g. "DISP"
    pub field_type: FieldType,
    /// Maps Node ID -> Vector of values (e.g. [dx, dy, dz])
    pub data: HashMap<NodeId, Vec<f64>>,
}

impl NodalResult {
    pub fn new(name: &str, field_type: FieldType) -> Self {
        Self {
            name: name.to_string(),
            field_type,
            data: HashMap::new(),
        }
    }

    pub fn insert(&mut self, node_id: NodeId, values: Vec<f64>) {
        self.data.insert(node_id, values);
    }
}

/// Holds all results for a single step (e.g. "Step 1 - Static")
#[allow(unused)]
#[derive(Debug, Clone)]
pub struct StepResult {
    pub step_id: usize,
    pub step_name: String,
    pub time_increment: f64,
    pub nodal_results: Vec<NodalResult>,
}

impl StepResult {
    pub fn new(step_id: usize, name: &str, time: f64) -> Self {
        Self {
            step_id,
            step_name: name.to_string(),
            time_increment: time,
            nodal_results: Vec::new(),
        }
    }
}
