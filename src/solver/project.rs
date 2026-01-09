use std::{
    collections::HashMap,
    error::Error,
    fmt::Write,
    fs::{self, read_to_string},
    path::PathBuf,
};

use crate::solver::{
    inp::InpFile,
    material::Material,
    mesh_lib::mesh::Mesh,
    parser::Parser,
    step::steps::Step,
};

/// This struct holds all relevant information from the `.inp` file
#[derive(Debug, Clone, Default)]
pub struct Project {
    /// The filepath to the .inp file.
    pub filepath: PathBuf,
    pub mesh: Box<Mesh>,
    pub steps: Vec<Step>,
    pub input: Box<InpFile>,
    pub materials: Vec<Material>,
    /// Map: Element-ID -> Material-Index
    pub element_materials: HashMap<usize, usize>,
}

impl Project {
    pub fn new() -> Self {
        Default::default()
    }

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

        if !path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("inp"))
        {
            let mut s = path.into_os_string();
            s.push(".inp");
            path = PathBuf::from(s);
        }

        let input_content = read_to_string(&path)?;
        if let Some(out) = preprocess_output {
            fs::write(out, &input_content)?;
        }
        
        let input = InpFile::new(&input_content);
        let mut project = Parser::new(&input).parse()?;
        
        project.filepath = path;
        project.input = Box::new(input);
        
        
        Ok(project)
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
