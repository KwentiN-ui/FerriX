//! Load and boundary condition definitions.
//!
//! This module defines concentrated loads and displacement boundary conditions
//! that can be applied to nodes within analysis steps.

use crate::solver::{amplitude::Amplitude, ids::NodeId, time::SolverTime};

/// A concentrated nodal load (*CLOAD).
#[derive(Debug, Clone)]
pub struct Load {
    node_id: NodeId,
    dof: usize, // 0=x, 1=y, 2=z
    value: f64,
    amplitude: Amplitude,
    origin_step: usize,
}

impl Load {
    /// Creates a new `Load`.
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

    /// Returns the ID of the node the load is applied to.
    #[must_use]
    pub fn node_id(&self) -> NodeId {
        self.node_id
    }

    /// Returns the degree of freedom (0, 1, or 2) the load is applied to.
    #[must_use]
    pub fn dof(&self) -> usize {
        self.dof
    }

    /// Calculates the instantaneous value of the load at a given simulation time.
    #[must_use]
    pub fn value(&self, time: &SolverTime, current_step: usize) -> f64 {
        self.amplitude.y(time, self.origin_step, current_step) * self.value
    }
}

/// A displacement boundary condition (*BOUNDARY).
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    node_id: NodeId,
    dof: usize, // 0=x, 1=y, 2=z
    value: f64,
    amplitude: Amplitude,
    origin_step: usize,
}

impl BoundaryCondition {
    /// Creates a new `BoundaryCondition`.
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

    /// Returns the ID of the node the boundary condition is applied to.
    #[must_use]
    pub fn node_id(&self) -> NodeId {
        self.node_id
    }

    /// Returns the degree of freedom (0, 1, or 2) being constrained.
    #[must_use]
    pub fn dof(&self) -> usize {
        self.dof
    }

    /// Calculates the instantaneous target value for the boundary condition.
    #[must_use]
    pub fn value(&self, time: &SolverTime, current_step: usize) -> f64 {
        self.amplitude.y(time, self.origin_step, current_step) * self.value
    }
}
