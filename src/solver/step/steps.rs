use crate::solver::error::Result;
use crate::solver::io::writer::ResultWriter;
use crate::solver::project::Project;
use crate::solver::state::SolutionState;
use crate::solver::step::static_step::StaticStep;
use crate::solver::time::SolverTime;
use strum::{EnumDiscriminants, EnumIter};

#[derive(Debug, Clone, EnumDiscriminants)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(StepKind))]
pub enum Step {
    StaticStep(StaticStep),
}

impl Step {
    /// Solves the current step.
    ///
    /// # Errors
    /// Returns an error if the step fails to solve.
    pub fn solve(
        &self,
        step_id: usize,
        project: &Project,
        solution_state: &mut SolutionState,
        writer: &dyn ResultWriter,
        timer: &mut SolverTime,
    ) -> Result<()> {
        match self {
            Step::StaticStep(static_step) => {
                timer.new_step(static_step.increment_data.time_period);
                static_step.solve(step_id, project, solution_state, writer, timer)
            }
        }
    }
}
