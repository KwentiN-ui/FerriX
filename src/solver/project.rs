use std::{
    error::Error,
    fmt::Write,
    fs::{self, read_to_string},
    path::PathBuf,
};

use crate::solver::{
    inp::InpFile,
    material::Material,
    mesh_lib::mesh::Mesh,
    section::{Section, SectionString},
    step::steps::Step,
};

/// This struct holds all relevant information from the `.inp` file
#[derive(Debug, Clone)]
pub struct Project {
    /// The filepath to the .inp file.
    pub filepath: PathBuf,
    pub mesh: Box<Mesh>,
    pub steps: Vec<Step>,
    pub input: Box<InpFile>,
    pub materials: Vec<Material>,
    pub sections: Vec<Section>,
}

impl Project {
    pub fn get_info(&self) -> String {
        let mut info = String::new();

        let _ = writeln!(info, "--- Project Info ---");
        let _ = writeln!(info, "Jobname:");
        let _ = writeln!(info, "  {}", self.jobname());
        let _ = writeln!(info, "Mesh:");
        let _ = writeln!(info, "  Nodes: {}", self.mesh.nodes.len());
        let _ = writeln!(info, "  Elements:");

        for (elem, count) in self.mesh.count_by_type() {
            let _ = writeln!(info, "  - {elem:?}: {count}");
        }

        info
    }

    pub fn from_jobname(
        jobname_filepath: &str,
        preprocess_output: Option<&String>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut path = PathBuf::from(jobname_filepath);

        // Appends .inp if not existing already
        if !path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("inp"))
        {
            let mut s = path.into_os_string();
            s.push(".inp");
            path = PathBuf::from(s);
        }

        // read and process the input
        let input: InpFile = InpFile::new(&read_to_string(&path)?);
        if let Some(out) = preprocess_output {
            fs::write(out, &input.0)?;
        }
        let inp_sections = InpSection::parse_sections_from_input(&input);

        let sections = dbg!(SectionString::parse_input(&input));

        let mesh = Mesh::from_sections(&input, &inp_sections)?;
        let materials = Material::from_input(&input);

        let steps = Step::parse_steps(&input);

        Ok(Self {
            steps,
            materials,
            sections,
            filepath: path,
            mesh: mesh.into(),
            input: input.into(),
        })
    }

    /// Gets the jobname from a filepath.
    pub fn jobname(&self) -> String {
        self.filepath
            .file_stem()
            .expect("The jobname could not be inferred by the given arguments.")
            .to_str()
            .expect("There is no reason why this should fail")
            .to_string()
    }
}

/// Different types of sections of the .inp file and their line-number
#[derive(Debug, PartialEq)]
pub enum InpSection {
    Node(usize),
    Element(usize),
    /// Stores: `LineNumber`, Name, `is_generate`
    Nset(usize, String, bool),
}

impl InpSection {
    fn parse_sections_from_input(file: &InpFile) -> Vec<InpSection> {
        let mut sections: Vec<InpSection> = Vec::new();
        for (nr, line) in file.0.lines().enumerate() {
            if line == "*NODE" {
                sections.push(InpSection::Node(nr));
            } else if line.starts_with("*ELEMENT") {
                sections.push(InpSection::Element(nr));
            } else if line.starts_with("*NSET") {
                let name = line
                    .split(',')
                    .find(|part| part.trim().starts_with("NSET="))
                    .map_or("UNKNOWN".to_string(), |part| {
                        part.split('=')
                            .nth(1)
                            .unwrap_or("UNKNOWN")
                            .trim()
                            .to_string()
                    });

                let is_generate = line
                    .chars()
                    .filter(|c| *c != ' ')
                    .collect::<String>()
                    .contains(",GENERATE");

                sections.push(InpSection::Nset(nr, name, is_generate));
            }
        }
        sections
    }
}
