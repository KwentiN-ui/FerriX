//! Time incrementation logic for simulation steps.
//!
//! This module defines how simulation time is advanced, including constraints
//! on minimum and maximum increment sizes.

/// Configuration for time incrementation in a simulation step.
#[derive(Debug, Clone, Copy)]
pub struct IncrementData {
    /// First time increment. Subsequent increments may be adjusted automatically.
    pub initial_time_increment: f64,
    /// Total duration of the simulation step.
    pub time_period: f64,
    /// Minimum allowed time increment. The simulation will fail if a smaller increment is required.
    pub min_time_increment: f64,
    /// Maximum allowed time increment. Useful for ensuring sufficient temporal resolution.
    pub max_time_increment: f64,
    /// Maximum number of increments allowed before the simulation is aborted.
    pub max_iterations: usize,
}

impl Default for IncrementData {
    fn default() -> Self {
        Self {
            initial_time_increment: 1.,
            time_period: 1.,
            min_time_increment: 1e-5,
            max_time_increment: 1.,
            max_iterations: 10_000,
        }
    }
}
