use crate::solver::ids::{NodeId, LoadId, BoundaryConditionId};

/// Represents a concentrated load (*CLOAD)
#[derive(Debug, Clone)]
pub struct Load {
    #[allow(dead_code)]
    pub id: LoadId,
    pub node_id: NodeId,
    pub dof: usize, // 0=x, 1=y, 2=z
    pub value: f64,
}

/// Represents a boundary condition (*BOUNDARY)
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    #[allow(dead_code)]
    pub id: BoundaryConditionId,
    pub node_id: NodeId,
    pub dof: usize, // 0=x, 1=y, 2=z
    pub value: f64,
}
