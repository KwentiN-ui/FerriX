//! Error handling for the `FerriX` solver.
//!
//! This module defines the custom error types used throughout the library to handle
//! IO, parsing, numerical, and convergence issues.

use crate::solver::ids::{ElementId, NodeId};
use std::path::PathBuf;
use thiserror::Error;

/// Custom error type for all `FerriX` operations.
#[derive(Error, Debug)]
pub enum FerrixError {
    /// Errors occurring during file input/output at a specific path.
    #[error("IO error at {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },

    /// General IO errors not tied to a specific path.
    #[error("Generic IO error: {0}")]
    GenericIo(#[from] std::io::Error),

    /// Errors encountered while parsing input files (e.g., .inp files).
    #[error("Parsing error at line {line}: {message}")]
    ParseError { line: usize, message: String },

    /// Errors for valid syntax keywords that are not yet implemented in `FerriX`.
    #[error("Unsupported keyword at line {line}: {keyword}")]
    UnsupportedKeyword { line: usize, keyword: String },

    /// Errors when an unrecognized solver algorithm is requested.
    #[error("Unknown solver type: {0}")]
    UnknownSolver(String),

    /// Error when a requested node is missing from the mesh.
    #[error("Node {0} not found in mesh")]
    NodeNotFound(NodeId),

    /// Error when a requested element is missing from the mesh.
    #[error("Element {0} not found in mesh")]
    ElementNotFound(ElementId),

    /// Error when a referenced material property is undefined.
    #[error("Material {0} not found")]
    MaterialNotFound(String),

    /// Error when a referenced element set is missing.
    #[error("Elset {0} not found")]
    ElsetNotFound(String),

    /// Error when a referenced node set is missing.
    #[error("Nset {0} not found")]
    NsetNotFound(String),

    /// Errors related to linear algebra or numerical stability.
    #[error("Mathematical error: {0}")]
    NumericalError(String),

    /// Error when iterative solvers fail to meet convergence criteria.
    #[error("Solver failed to converge: {0}")]
    ConvergenceError(String),

    /// Error for invalid model configurations or inconsistent data.
    #[error("Invalid model state: {0}")]
    InvalidModelState(String),

    /// General catch-all error with a message.
    #[error("Generic error: {0}")]
    Generic(String),
}

/// A specialized Result type for `FerriX` operations.
pub type Result<T> = std::result::Result<T, FerrixError>;
