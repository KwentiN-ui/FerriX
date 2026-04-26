use crate::solver::{amplitude::Amplitude, ids::NodeId, time::SolverTime};

/// Represents a concentrated load (*CLOAD)
#[derive(Debug, Clone)]
pub struct Load {
    node_id: NodeId,
    dof: usize, // 0=x, 1=y, 2=z
    value: f64,
    amplitude: Amplitude,
    origin_step: usize,
}

impl Load {
    #[must_use] 
    pub fn new(
        node_id: NodeId,
        dof: usize,
        value: f64,
        amplitude: Option<Amplitude>,
        origin_step: usize,
    ) -> Self {
        Self {
            node_id,
            dof,
            value,
            amplitude: amplitude.unwrap_or_default(),
            origin_step,
        }
    }
    #[must_use] 
    pub fn node_id(&self) -> NodeId {
        self.node_id
    }
    #[must_use] 
    pub fn dof(&self) -> usize {
        self.dof
    }
    #[must_use] 
    pub fn value(&self, time: &SolverTime, current_step: usize) -> f64 {
        self.amplitude.y(time, self.origin_step, current_step) * self.value
    }
}

/// Represents a boundary condition (*BOUNDARY)
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    node_id: NodeId,
    dof: usize, // 0=x, 1=y, 2=z
    value: f64,
    amplitude: Amplitude,
    origin_step: usize,
}

impl BoundaryCondition {
    #[must_use] 
    pub fn new(
        node_id: NodeId,
        dof: usize,
        value: f64,
        amplitude: Option<Amplitude>,
        origin_step: usize,
    ) -> Self {
        Self {
            node_id,
            dof,
            value,
            amplitude: amplitude.unwrap_or_default(),
            origin_step,
        }
    }
    #[must_use] 
    pub fn node_id(&self) -> NodeId {
        self.node_id
    }
    #[must_use] 
    pub fn dof(&self) -> usize {
        self.dof
    }
    #[must_use] 
    pub fn value(&self, time: &SolverTime, current_step: usize) -> f64 {
        self.amplitude.y(time, self.origin_step, current_step) * self.value
    }
}
