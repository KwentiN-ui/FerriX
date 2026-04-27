//! Global simulation state management.
//!
//! This module defines the `SolutionState` struct, which holds the cumulative
//! results (like displacements) as the simulation progresses through steps.

/// Represents the current physical state of the entire model.
///
/// This struct is passed between analysis steps, carrying the cumulative
/// effects of previous computations and providing a starting point for the next.
#[derive(Debug, Clone)]
pub struct SolutionState {
    /// Cumulative nodal displacements at the end of the last completed increment.
    /// This vector is indexed by the global degree of freedom (3 per node).
    pub displacements: Vec<f64>,
    /// Nodal temperatures. Vector is indexed by node index.
    pub temperatures: Vec<f64>,
}

impl SolutionState {
    /// Creates a new, empty `SolutionState` for the beginning of an analysis.
    #[must_use]
    pub fn new(num_dofs: usize, num_nodes: usize) -> Self {
        Self {
            displacements: vec![0.0; num_dofs],
            temperatures: vec![0.0; num_nodes],
        }
    }
}
