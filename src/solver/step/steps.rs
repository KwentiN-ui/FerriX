use strum::{EnumDiscriminants, EnumIter};

use crate::solver::io::writer::ResultWriter;
use crate::solver::project::Project;
use crate::solver::state::SolutionState;
use crate::solver::step::static_step::StaticStep;
use std::error::Error;

#[derive(Debug, Clone, EnumDiscriminants)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(StepKind))]
pub enum Step {
    StaticStep(StaticStep),
}

impl Step {
    pub fn solve(
        &self,
        step_id: usize,
        project: &Project,
        solution_state: &mut SolutionState,
        writer: &mut dyn ResultWriter,
    ) -> Result<(), Box<dyn Error>> {
        match self {
            Step::StaticStep(static_step) => {
                static_step.solve(step_id, project, solution_state, writer)
            }
        }
    }
}
