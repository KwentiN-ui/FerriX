//! Abaqus-style input file (.inp) handling.
//!
//! This module provides structures for reading and preprocessing input files,
//! stripping comments and normalizing whitespace before parsing.

use derive_more::{AsRef, Deref, Display};

use crate::solver::parser::preprocess_inp;

/// Represents a preprocessed input file.
#[derive(Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Display, Deref, AsRef)]
pub struct InpFile(pub String);

impl InpFile {
    /// Creates a new `InpFile` by preprocessing the raw input string.
    #[must_use]
    pub fn new(input: &str) -> Self {
        InpFile(preprocess_inp(input))
    }
}
