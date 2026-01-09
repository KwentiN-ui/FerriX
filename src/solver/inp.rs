use crate::solver::parsing::preprocess_inp;

#[derive(Debug, Clone)]
pub struct InpFile(pub String);

impl InpFile {
    pub fn new(input: &str) -> Self {
        InpFile(preprocess_inp(input))
    }
}
