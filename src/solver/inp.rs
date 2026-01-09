use crate::solver::parser::preprocess_inp;

#[derive(Debug, Clone, Default)]
pub struct InpFile(pub String);

impl InpFile {
    pub fn new(input: &str) -> Self {
        InpFile(preprocess_inp(input))
    }
}
