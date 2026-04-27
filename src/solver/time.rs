//! Simulation time tracking.
//!
//! This module provides the `SolverTime` utility, which maintains both the
//! global simulation time and the local time relative to the current step.

/// Tracks global and local simulation time across steps and increments.
#[derive(Debug, Clone)]
pub struct SolverTime {
    /// Cumulative simulation time from the beginning of the project.
    global: f64,
    /// Time elapsed within the current analysis step.
    local: f64,
    /// Total duration of the current analysis step.
    local_max: f64,
}

impl Default for SolverTime {
    fn default() -> Self {
        Self::new()
    }
}

impl SolverTime {
    /// Creates a new `SolverTime` initialized to zero.
    #[must_use]
    pub fn new() -> Self {
        Self {
            global: 0.,
            local: 0.,
            local_max: 1.,
        }
    }

    /// Resets the local step-time for a new simulation step.
    pub fn new_step(&mut self, max_time: f64) {
        self.local = 0.;
        self.local_max = max_time;
    }

    /// Advances both global and local time by the given increment.
    pub fn new_increment(&mut self, timestep: f64) {
        self.global += timestep;
        self.local += timestep;
    }

    /// Returns the current local time relative to the start of the step.
    #[must_use]
    pub fn local_time(&self) -> f64 {
        self.local
    }

    /// Returns the cumulative global simulation time.
    #[must_use]
    pub fn global_time(&self) -> f64 {
        self.global
    }

    /// Returns the maximum duration of the current step.
    #[must_use]
    pub fn local_max_time(&self) -> f64 {
        self.local_max
    }
}
