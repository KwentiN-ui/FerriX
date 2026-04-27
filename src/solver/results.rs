//! Simulation results and field data.
//!
//! This module defines the structures used to store and organize results from
//! different increments and steps, such as displacements, stresses, and strains.

use crate::solver::ids::NodeId;
use std::collections::HashMap;

/// Types of physical fields that can be stored in results.
#[allow(unused)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldType {
    /// Nodal displacements (U).
    Displacement,
    /// Element stresses (S).
    Stress,
    /// Element strains (E).
    Strain,
}

/// Stores result values for a specific field across all relevant nodes.
#[derive(Debug, Clone)]
pub struct NodalResult {
    /// Descriptive name of the result field (e.g., "DISP").
    #[allow(dead_code)]
    pub name: String,
    /// The type of field being stored.
    pub field_type: FieldType,
    /// Mapping from `NodeId` to the vector of field values (e.g., [dx, dy, dz]).
    pub data: HashMap<NodeId, Vec<f64>>,
}

impl NodalResult {
    /// Creates a new, empty `NodalResult` container.
    #[must_use]
    pub fn new(name: &str, field_type: FieldType) -> Self {
        Self {
            name: name.to_string(),
            field_type,
            data: HashMap::new(),
        }
    }

    /// Inserts field values for a specific node.
    pub fn insert(&mut self, node_id: NodeId, values: Vec<f64>) {
        self.data.insert(node_id, values);
    }
}

/// Aggregates all results for a single simulation increment.
#[allow(unused)]
#[derive(Debug, Clone)]
pub struct IncResult {
    /// Index of the current step.
    pub step_id: usize,
    /// Index of the current increment within the step.
    pub inc_id: usize,
    /// Name of the current step.
    pub step_name: String,
    /// Time increment size for this state.
    pub time_increment: f64,
    /// Collection of nodal results for this increment.
    pub nodal_results: Vec<NodalResult>,
}

impl IncResult {
    #[must_use]
    pub fn new(step_id: usize, inc_id: usize, name: &str, time: f64) -> Self {
        Self {
            step_id,
            inc_id,
            step_name: name.to_string(),
            time_increment: time,
            nodal_results: Vec::new(),
        }
    }
}
