//! This module contains the `SolutionState` struct, which holds the cumulative
//! state of the simulation at a given point in time.

/// Holds the cumulative state of the simulation.
///
/// This struct is passed from one analysis step to the next, carrying the results
/// and allowing steps to build upon each other.
#[derive(Debug, Clone)]
pub struct SolutionState {
    /// Cumulative nodal displacements at the end of the last completed increment.
    /// This vector is indexed by the global degree of freedom.
    pub displacements: Vec<f64>,
}

impl SolutionState {
    /// Creates a new, empty `SolutionState` for the beginning of an analysis.
    #[must_use] 
    pub fn new(num_dofs: usize) -> Self {
        Self {
            displacements: vec![0.0; num_dofs],
        }
    }
}
