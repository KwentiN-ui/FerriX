/// Keeps track of the current simulation time.
/// The program differentiates into "local" step-time and global "simulation time"
#[derive(Debug, Clone)]
pub struct SolverTime {
    global: f64,
    local: f64,
    /// The timespan of the currently running step
    local_max: f64,
}

impl SolverTime {
    pub fn new() -> Self {
        Self {
            global: 0.,
            local: 0.,
            local_max: 1.,
        }
    }
    pub fn new_step(&mut self, max_time: f64) {
        self.local = 0.;
        self.local_max = max_time;
    }
    pub fn new_increment(&mut self, timestep: f64) {
        self.global += timestep;
        self.local += timestep;
    }
    pub fn local_time(&self) -> f64 {
        self.local
    }
    pub fn global_time(&self) -> f64 {
        self.global
    }
    pub fn local_max_time(&self) -> f64 {
        self.local_max
    }
}
