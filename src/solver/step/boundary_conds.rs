use crate::solver::{amplitude::Amplitude, ids::NodeId, time::SolverTime};

/// Represents a concentrated load (*CLOAD)
#[derive(Debug, Clone)]
pub struct Load {
    node_id: NodeId,
    dof: usize, // 0=x, 1=y, 2=z
    value: f64,
    amplitude: Amplitude,
}

impl Load {
    pub fn new(node_id: NodeId, dof: usize, value: f64, amplitude: Option<Amplitude>) -> Self {
        Self {
            node_id,
            dof,
            value,
            amplitude: amplitude.unwrap_or_default(),
        }
    }
    pub fn node_id(&self) -> NodeId {
        self.node_id
    }
    pub fn dof(&self) -> usize {
        self.dof
    }
    pub fn value(&self, time: &SolverTime) -> f64 {
        self.amplitude.y(time) * self.value
    }
}

/// Represents a boundary condition (*BOUNDARY)
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    node_id: NodeId,
    dof: usize, // 0=x, 1=y, 2=z
    value: f64,
    amplitude: Amplitude,
}

impl BoundaryCondition {
    pub fn new(node_id: NodeId, dof: usize, value: f64, amplitude: Option<Amplitude>) -> Self {
        Self {
            node_id,
            dof,
            value,
            amplitude: amplitude.unwrap_or_default(),
        }
    }
    pub fn node_id(&self) -> NodeId {
        self.node_id
    }
    pub fn dof(&self) -> usize {
        self.dof
    }
    pub fn value(&self, time: &SolverTime) -> f64 {
        self.amplitude.y(time) * self.value
    }
}
