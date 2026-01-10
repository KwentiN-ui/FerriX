//! This module contains all custom ID types for better type safety at compile-time.

use derive_more::{Deref, Display, From, Into};

/// Node ID as defined in the INP File
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct NodeId(pub usize);

/// Element ID as defined in the INP File
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct ElementId(pub usize);

/// Unique identifier for a `Load`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct LoadId(pub usize);

/// Unique identifier for a `BoundaryCondition`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct BoundaryConditionId(pub usize);
