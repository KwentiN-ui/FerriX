use crate::solver::ids::{ElementId, NodeId};
use std::path::PathBuf;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FerrixError {
    #[error("IO error at {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Generic IO error: {0}")]
    GenericIo(#[from] std::io::Error),

    #[error("Parsing error at line {line}: {message}")]
    ParseError { line: usize, message: String },

    #[error("Unsupported keyword at line {line}: {keyword}")]
    UnsupportedKeyword { line: usize, keyword: String },

    #[error("Unknown solver type: {0}")]
    UnknownSolver(String),

    #[error("Node {0} not found in mesh")]
    NodeNotFound(NodeId),

    #[error("Element {0} not found in mesh")]
    ElementNotFound(ElementId),

    #[error("Material {0} not found")]
    MaterialNotFound(String),

    #[error("Elset {0} not found")]
    ElsetNotFound(String),

    #[error("Nset {0} not found")]
    NsetNotFound(String),

    #[error("Mathematical error: {0}")]
    NumericalError(String),

    #[error("Solver failed to converge: {0}")]
    ConvergenceError(String),

    #[error("Invalid model state: {0}")]
    InvalidModelState(String),

    #[error("Generic error: {0}")]
    Generic(String),
}

pub type Result<T> = std::result::Result<T, FerrixError>;
