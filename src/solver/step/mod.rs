//! This submodule houses step-related logic. Every step-type can fully define it's own parsing logic,
//! as well as the solving procedure. New steps need to be included in the `Step` enum. The `main` method then
//! calls subsequent steps that were found through TODO

pub mod boundary_conds;
pub mod static_step;
pub mod steps;
