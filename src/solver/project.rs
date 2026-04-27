//! FEA project management and central data coordination.
//!
//! The `Project` struct in this module acts as the "God Object" or central hub,
//! aggregating all mesh, material, and step data into a single container.

use crate::solver::error::{FerrixError, Result};
use crate::solver::{amplitude::Amplitude, ids::ElementId};
use std::{
    collections::HashMap,
    fmt::Write,
    fs::{self, read_to_string},
    path::PathBuf,
};

use crate::solver::inp::InpFile;
use crate::solver::material::Material;
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::parser::Parser;
use crate::solver::step::{
    boundary_conds::{BoundaryCondition, Load},
    steps::Step,
};
use std::sync::Arc;

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
    pub materials: Vec<Arc<dyn Material>>,
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
        let _ = writeln!(
            info,
            "  {}",
            self.jobname().unwrap_or_else(|_| "Unknown".to_string())
        );
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
    ) -> Result<Self> {
        let mut path = PathBuf::from(jobname_filepath);

        if !path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("inp"))
        {
            let mut s = path.into_os_string();
            s.push(".inp");
            path = PathBuf::from(s);
        }

        let input_content = read_to_string(&path).map_err(|e| FerrixError::Io {
            path: path.clone(),
            source: e,
        })?;
        if let Some(out) = preprocess_output {
            let out_path = PathBuf::from(out);
            fs::write(&out_path, &input_content).map_err(|e| FerrixError::Io {
                path: out_path,
                source: e,
            })?;
        }

        let input = InpFile::new(&input_content);
        let mut project = Parser::new(&input).parse()?;

        project.filepath = path;
        project.input = Box::new(input);

        Ok(project)
    }

    /// Gets the jobname from a filepath.
    ///
    /// # Errors
    /// Returns an error if the jobname cannot be inferred.
    pub fn jobname(&self) -> Result<String> {
        let stem = self
            .filepath
            .file_stem()
            .ok_or_else(|| {
                FerrixError::InvalidModelState("Could not infer jobname from path".into())
            })?
            .to_str()
            .ok_or_else(|| FerrixError::InvalidModelState("Invalid UTF-8 in jobname".into()))?;
        Ok(stem.to_string())
    }

    /// Filepath to the output directory
    ///
    /// # Errors
    /// Returns an error if the job directory cannot be inferred.
    pub fn job_dir(&self) -> Result<PathBuf> {
        self.filepath
            .parent()
            .map(PathBuf::from)
            .ok_or_else(|| FerrixError::InvalidModelState("Could not infer job directory".into()))
    }
}
