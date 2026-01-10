use strum::{EnumDiscriminants, EnumIter};

use crate::solver::solvers::SolverType;

#[derive(Debug, Clone, EnumDiscriminants)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(StepKind))]
pub enum Step {
    StaticStep(SolverType),
}

impl Step {}
