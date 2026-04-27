//! Type-safe identifiers for FEA entities.
//!
//! This module provides wrapper types for IDs (nodes, elements, loads, etc.) to prevent
//! accidental mixing of different identifier types at compile-time.

use derive_more::{Deref, Display, From, Into};

/// A unique identifier for a node in the mesh, typically mapping to the ID in the input file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct NodeId(pub usize);

/// A unique identifier for an element in the mesh, typically mapping to the ID in the input file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct ElementId(pub usize);

/// A unique identifier for a load applied to the model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct LoadId(pub usize);

/// A unique identifier for a boundary condition applied to the model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, From, Into, Deref)]
pub struct BoundaryConditionId(pub usize);
