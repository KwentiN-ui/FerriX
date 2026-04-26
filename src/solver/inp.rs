use derive_more::{AsRef, Deref, Display};

use crate::solver::parser::preprocess_inp;

#[derive(Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Display, Deref, AsRef)]
pub struct InpFile(pub String);

impl InpFile {
    #[must_use] 
    pub fn new(input: &str) -> Self {
        InpFile(preprocess_inp(input))
    }
}
