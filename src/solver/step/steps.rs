use strum_macros::EnumIter;

/// Each possible Step type is listed here. They contain the corresponding linenumber to ensure the right step is read.
#[derive(Debug, Clone, EnumIter)]
pub enum Step {
    StaticStep(usize),
}

impl Step {
    /// If the step is to be included in the analysis, then
    /// at least one line in the input file needs to pass `line.starts_with(step.keyword())`.
    pub fn keyword(&self) -> &'static str {
        match self {
            Step::StaticStep(_) => "*STATIC",
        }
    }
}
