#[derive(Debug, Clone, Copy)]
pub struct IncrementData {
    /// First time increment. Subsequent increments are chosen automatically
    pub initial_time_increment: f64,
    /// Length of the timestep
    pub time_period: f64,
    /// Minimum allowed time increment
    pub min_time_increment: f64,
    /// Maximum allowed time increment. Choose a lower amount if you want more result files
    pub max_time_increment: f64,
    /// The simulation aborts if this number of increments is reached
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
