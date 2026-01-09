use strum::{EnumDiscriminants, EnumIter, IntoEnumIterator};

use crate::solver::inp::InpFile;

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

impl Step {
    /// This functions searches the input file for known step-types. These are later run in the order defined here.
    pub fn parse_steps(input: &InpFile) -> Vec<Step> {
        let mut steps = Vec::new();
        let mut lines = input.0.lines().enumerate().peekable();

        while let Some((_, line)) = lines.next() {
            if line.starts_with("*STEP") {
                if let Some((next_nr, next_line)) = lines.peek() {
                    for step in StepKind::iter() {
                        if next_line.starts_with(step.keyword()) {
                            steps.push(step.create(*next_nr));
                            break;
                        }
                    }
                }
            }
        }
        steps
    }
}
