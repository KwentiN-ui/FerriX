use strum::{EnumDiscriminants, EnumIter, IntoEnumIterator};

use crate::solver::inp::InpFile;

#[derive(Debug, Clone, EnumDiscriminants, PartialEq)]
#[strum_discriminants(derive(EnumIter))]
#[strum_discriminants(name(SectionString))]
pub enum Section {
    /// Solid Section Definition ()
    SolidSection(String, String),
}

impl SectionString {
    pub fn parse_input(input: &InpFile) -> Vec<Section> {
        let mut sections = Vec::new();

        let mut lines = input.0.lines().peekable();
        while let Some(line) = lines.next() {
            for sec in SectionString::iter() {
                if line.starts_with(sec.keyword()) {
                    sections.push(sec.create(line, lines.peek().map(|v| &**v)));
                }
            }
        }
        sections
    }

    pub fn keyword(&self) -> &str {
        match self {
            SectionString::SolidSection => "*SOLID SECTION",
        }
    }

    pub fn create(self, line: &str, _next_line: Option<&str>) -> Section {
        match self {
            SectionString::SolidSection => {
                let args: Vec<&str> = line.split(',').skip(1).map(str::trim).collect();

                Section::SolidSection(
                    args[0].split('=').next_back().unwrap().trim().to_string(),
                    args[1].split('=').next_back().unwrap().trim().to_string(),
                )
            }
        }
    }
}

#[cfg(test)]
mod test {
    use itertools::assert_equal;

    use super::*;

    #[test]
    fn test_input_parsing() {
        let input = InpFile::new(
            "*Solid section, Elset=Internal_Selection-1_Solid_Section-1, Material=S235
            *Solid section, Elset=Internal_Selection-1_Solid_Section-2, Material=S355",
        );

        let parsed = SectionString::parse_input(&input);

        let correct = vec![
            Section::SolidSection(
                "INTERNAL_SELECTION-1_SOLID_SECTION-1".to_string(),
                "S235".to_string(),
            ),
            Section::SolidSection(
                "INTERNAL_SELECTION-1_SOLID_SECTION-2".to_string(),
                "S355".to_string(),
            ),
        ];

        assert_equal(correct, parsed);
    }
}
