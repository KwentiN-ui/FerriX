/// Represents a concentrated load (*CLOAD)
#[derive(Debug, Clone)]
pub struct Load {
    pub node_id: usize,
    pub dof: usize, // 0=x, 1=y, 2=z
    pub value: f64,
}

/// Represents a boundary condition (*BOUNDARY)
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    pub node_id: usize,
    pub dof: usize, // 0=x, 1=y, 2=z
    pub value: f64,
}
