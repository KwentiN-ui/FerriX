//! This module contains all custom ID types for better type safety at compile-time.

use std::fmt;

// TODO differentiate between local and global IDs
/// A typesafe Node ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct NodeId(pub usize);

impl fmt::Display for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A typesafe Element ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ElementId(pub usize);

impl fmt::Display for ElementId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A typesafe Load ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LoadId(pub usize);

impl fmt::Display for LoadId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A typesafe `BoundaryCondition` ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryConditionId(pub usize);

impl fmt::Display for BoundaryConditionId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
