use strum::{EnumDiscriminants, EnumIter};

#[derive(Debug, Clone, EnumDiscriminants)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(StepKind))]
pub enum Step {
    StaticStep(usize),
}

impl StepKind {
    /// If a `*STEP` card is followed by this keyword, it will be included in the analysis.
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
