use strum::{EnumDiscriminants, EnumIter};

#[derive(Debug, Clone, EnumDiscriminants)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(StepKind))]
pub enum Step {
    StaticStep(usize),
}

impl StepKind {
    /// If the step is to be included in the analysis, then
    /// at least one line in the input file needs to pass `line.starts_with(step.keyword())`.
    pub fn keyword(self) -> &'static str {
        match self {
            StepKind::StaticStep => "*STATIC",
        }
    }
    pub fn create(self, line_number: usize) -> Step {
        match self {
            StepKind::StaticStep => Step::StaticStep(line_number),
        }
    }
}
