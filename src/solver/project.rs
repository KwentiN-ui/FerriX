use crate::solver::{amplitude::Amplitude, ids::ElementId};
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
    step::{
        boundary_conds::{BoundaryCondition, Load},
        steps::Step,
    },
};

/// Represents the central data structure for a single FEA problem.
///
/// This struct aggregates all the necessary information parsed from an `.inp` file,
/// including the mesh, materials, steps, loads, and boundary conditions. It serves
/// as the primary container passed between different components of the solver.
#[derive(Debug, Clone, Default)]
pub struct Project {
    /// The filepath to the .inp file, which serves as the job's identifier.
    pub filepath: PathBuf,
    /// The finite element mesh, containing all nodes, elements, and sets.
    pub mesh: Box<Mesh>,
    /// A vector of analysis steps to be executed in sequence.
    pub steps: Vec<Step>,
    /// The preprocessed content of the input file, stored for reference.
    pub input: Box<InpFile>,
    /// A list of all materials defined in the model.
    pub materials: Vec<Material>,
    /// A map that links each `ElementId` to its corresponding material's index in the `materials` vector.
    pub element_materials: HashMap<ElementId, usize>,
    /// A collection of all concentrated loads (*CLOAD) defined before the first step.
    pub initial_loads: Vec<Load>,
    /// A collection of all boundary conditions (*BOUNDARY) defined before the first step.
    pub initial_bcs: Vec<BoundaryCondition>,
    /// Nodal output variables requested by the user (e.g. `U`, `RF`).
    pub nodal_output: Vec<String>,
    /// Element output variables requested by the user (e.g. `S`, `E`).
    pub element_output: Vec<String>,
    /// Amplitude definitions as defined in the input file.
    pub amplitudes: HashMap<String, Amplitude>,
}

impl Project {
    #[must_use]
    pub fn new() -> Self {
        Project::default()
    }

    #[must_use]
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

    /// Creates a new Project from a jobname.
    ///
    /// # Errors
    /// Returns an error if the input file cannot be read or parsed.
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
    ///
    /// # Panics
    /// Panics if the jobname cannot be inferred.
    #[must_use]
    pub fn jobname(&self) -> String {
        self.filepath
            .file_stem()
            .expect("The jobname could not be inferred by the given arguments.")
            .to_str()
            .expect("There is no reason why this should fail")
            .to_string()
    }

    /// Filepath to the output directory
    ///
    /// # Panics
    /// Panics if the job directory cannot be inferred.
    #[must_use]
    pub fn job_dir(&self) -> PathBuf {
        self.filepath.parent().expect("Invalid filepath.").into()
    }
}
