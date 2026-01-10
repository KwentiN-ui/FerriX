//! This module contains all custom ID types for better type safety at compile-time.

use derive_more::{Deref, Display, From, Into};

// TODO differentiate between local and global IDs
/// A typesafe Node ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct NodeId(pub usize);

/// A typesafe Element ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct ElementId(pub usize);

/// A typesafe Load ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct LoadId(pub usize);

/// A typesafe `BoundaryCondition` ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct BoundaryConditionId(pub usize);
