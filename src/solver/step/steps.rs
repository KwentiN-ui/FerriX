use crate::solver::step::static_step::StaticStep;

#[derive(Debug, Clone)]
pub enum Step {
    StaticStep(StaticStep),
}
