//! Core solver module for Finite Element Analysis.
//!
//! This module organizes the various components required to solve FEA problems,
//! including mesh management, material properties, step definitions, and linear solvers.

pub mod amplitude;
pub mod assembler;
pub mod error;
pub mod ids;
pub mod increment;
pub mod inp;
pub mod io;
pub mod material;
pub mod mesh_lib;
pub mod parser;
pub mod preconditioner;
pub mod project;
pub mod results;
pub mod solvers;
pub mod state;
pub mod step;
pub mod time;
