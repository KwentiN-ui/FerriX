//! Global simulation state management.
//!
//! This module defines the `SolutionState` struct, which holds the cumulative
//! results (like displacements) as the simulation progresses through steps.

use crate::solver::ids::ElementId;
use crate::solver::project::Project;
use std::collections::HashMap;

/// Represents the current physical state of the entire model.
///
/// This struct is passed between analysis steps, carrying the cumulative
/// effects of previous computations and providing a starting point for the next.
#[derive(Debug, Clone, Default)]
pub struct SolutionState {
    /// Cumulative nodal displacements at the end of the last completed increment.
    /// This vector is indexed by the global degree of freedom (3 per node).
    pub displacements: Vec<f64>,
    /// Nodal temperatures at the start of the analysis (initial state).
    pub initial_temperatures: Vec<f64>,
    /// Current nodal temperatures.
    pub temperatures: Vec<f64>,
    /// State-dependent variables (SDVs) for each element at each integration point.
    /// Map: `ElementId` -> Vec<Vec<`SDV_Values`>> (outer Vec is IP, inner Vec is SDVs)
    pub material_states: HashMap<ElementId, Vec<Vec<f64>>>,
}

impl SolutionState {
    /// Creates a new, empty `SolutionState` for the beginning of an analysis.
    #[must_use]
    pub fn new(num_dofs: usize, num_nodes: usize) -> Self {
        Self {
            displacements: vec![0.0; num_dofs],
            initial_temperatures: vec![0.0; num_nodes],
            temperatures: vec![0.0; num_nodes],
            material_states: HashMap::new(),
        }
    }

    /// Initializes the solution state from a project.
    ///
    /// This sets up initial temperatures and material states (SDVs).
    ///
    /// # Panics
    /// Panics if an element in the mesh does not have a corresponding material assignment.
    pub fn initialize(&mut self, project: &Project) {
        // Initialize temperatures with default project-wide temperature
        self.initial_temperatures
            .fill(project.default_initial_temperature);
        self.temperatures.fill(project.default_initial_temperature);

        // Apply specific nodal initial temperatures from *INITIAL CONDITIONS
        for (&node_id, &temp) in &project.initial_temperatures {
            if let Some(idx) = project.mesh.get_index_for_node_id(node_id) {
                self.initial_temperatures[idx] = temp;
                self.temperatures[idx] = temp;
            }
        }

        // Initialize Material States (SDVs)
        for element in project.mesh.elements.values() {
            let material_index = project
                .element_materials
                .get(&element.get_id())
                .expect("Element has no material assigned");
            let material = &project.materials[*material_index];
            let num_sdvs = material.num_state_variables();
            let num_ips = element.integration_points().len();

            if num_sdvs > 0 {
                self.material_states
                    .insert(element.get_id(), vec![vec![0.0; num_sdvs]; num_ips]);
            }
        }
    }
}
