use strum::{EnumDiscriminants, EnumIter};

#[derive(Debug, Clone, EnumDiscriminants)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(StepKind))]
pub enum Step {
    StaticStep(usize),
}

impl Step {
}
