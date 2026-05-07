//! Input file parsing engine.
//!
//! This module implements the parser for Abaqus-style `.inp` files. It converts
//! the preprocessed text into a structured `Project` containing all FEA data.

use crate::solver::amplitude::{Amplitude, TimeSeries};
use crate::solver::error::{FerrixError, Result};
use crate::solver::ids::{ElementId, NodeId};
use crate::solver::increment::IncrementData;
use crate::solver::inp::InpFile;
use crate::solver::material::{BaseMaterial, TemperatureDependentLUT};
use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::mesh_lib::node::Node;
use crate::solver::project::Project;
use crate::solver::solvers::SolverType;
use crate::solver::step::boundary_conds::{BoundaryCondition, Load};
use crate::solver::step::static_step::StaticStep;
use crate::solver::step::steps::Step;
use std::collections::HashMap;
use std::str::FromStr;
use std::sync::Arc;
use strum_macros::EnumString;

/// Recognized keywords (cards) in the `.inp` file format.
#[derive(Debug, EnumString, PartialEq, Clone, Copy)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Keyword {
    Node,
    Element,
    Nset,
    Elset,
    Material,
    Elastic,
    Density,
    SolidSection,
    Step,
    Static,
    Boundary,
    Cload,
    EndStep,
    Heading,
    NodeFile,
    ElFile,
    Amplitude,
    Surface,
    Dload,
    Dsload,
    NodePrint,
    ElPrint,
    Include,
    BeamSection,
    Expansion,
    ShellSection,
    Conductivity,
    SpecificHeat,
    InitialConditions,
    Depvar,
}

/// The main parser for converting input file text into a `Project` model.
pub struct Parser<'a> {
    input: &'a InpFile,
    project: Project,
    materials: Vec<BaseMaterial>,
    current_keyword: Option<Keyword>,
    line_nr: usize,
    // Parser state
    element_type: Option<String>,
    elset_name: Option<String>,
    set_name: Option<String>,
    is_generate: bool,
    load_id_counter: usize,
    bc_id_counter: usize,
    current_step: Option<StaticStep>,
    step_counter: usize,
    step_loads: Vec<Load>,
    step_bcs: Vec<BoundaryCondition>,
    // keyword specific state
    initial_condition_type: Option<String>,
}

impl<'a> Parser<'a> {
    /// Creates a new `Parser` for the given input file.
    #[must_use]
    pub fn new(input: &'a InpFile) -> Self {
        Self {
            input,
            project: Project::new(),
            materials: Vec::new(),
            current_keyword: None,
            line_nr: 0,
            element_type: None,
            elset_name: None,
            set_name: None,
            is_generate: false,
            load_id_counter: 0,
            bc_id_counter: 0,
            current_step: None,
            step_counter: 0,
            step_loads: Vec::new(),
            step_bcs: Vec::new(),
            initial_condition_type: None,
        }
    }

    /// Parses the input file and returns a `Project`.
    ///
    /// # Errors
    /// Returns `FerrixError` if parsing fails due to syntax errors or unsupported features.
    pub fn parse(mut self) -> Result<Project> {
        let mut lines = self.input.lines().enumerate().peekable();
        while let Some((line_nr, line_content)) = lines.next() {
            self.line_nr = line_nr + 1;

            if line_content.starts_with('*') {
                self.parse_keyword(line_content, &mut lines)?;
                continue;
            }

            let line_content = line_content.trim();
            if line_content.is_empty() {
                continue;
            }

            if self.current_keyword.is_some() {
                self.parse_data(line_content)?;
            }
        }

        // Finalize last step if END STEP was missing
        if let Some(mut step) = self.current_step.take() {
            step.loads.clone_from(&self.step_loads);
            step.bcs.clone_from(&self.step_bcs);
            self.project.steps.push(Step::StaticStep(step));
        }

        // Move materials to project
        for mat in self.materials {
            self.project.materials.push(Arc::new(mat));
        }

        // Post-parsing steps, e.g. building node mappings
        self.project.mesh.build_node_mappings();

        Ok(self.project)
    }

    #[allow(clippy::too_many_lines)]
    fn parse_keyword(
        &mut self,
        line: &str,
        lines: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines>>,
    ) -> Result<()> {
        let parts: Vec<&str> = line.split(',').map(str::trim).collect();
        let keyword_str = parts[0]
            .strip_prefix('*')
            .unwrap_or("")
            .trim()
            .to_uppercase();

        if let Ok(keyword) = Keyword::from_str(&keyword_str.replace(' ', "_")) {
            self.current_keyword = Some(keyword);
            match keyword {
                Keyword::Element => {
                    let kwargs = get_keyword_arguments(line);
                    self.element_type = kwargs
                        .get("TYPE")
                        .and_then(Option::as_deref)
                        .map(str::to_string);
                    self.elset_name = kwargs
                        .get("ELSET")
                        .and_then(Option::as_deref)
                        .map(str::to_string);
                }
                Keyword::Nset => {
                    let kwargs = get_keyword_arguments(line);
                    self.set_name = kwargs
                        .get("NSET")
                        .and_then(Option::as_deref)
                        .map(str::to_string);
                    self.is_generate = kwargs.contains_key("GENERATE");

                    while let Some((_, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }
                        let (_, line_content) = lines.next().unwrap();
                        self.parse_nset(line_content);
                    }
                }
                Keyword::Elset => {
                    let kwargs = get_keyword_arguments(line);
                    self.set_name = kwargs
                        .get("ELSET")
                        .and_then(Option::as_deref)
                        .map(str::to_string);
                    self.is_generate = kwargs.contains_key("GENERATE");

                    while let Some((_, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }
                        let (_, line_content) = lines.next().unwrap();
                        self.parse_elset(line_content);
                    }
                }
                Keyword::Material => {
                    let kwargs = get_keyword_arguments(line);
                    let name = kwargs
                        .get("NAME")
                        .and_then(Option::as_deref)
                        .ok_or_else(|| FerrixError::ParseError {
                            line: self.line_nr,
                            message: "Material name not found".into(),
                        })?;
                    self.materials.push(BaseMaterial {
                        name: name.to_string(),
                        density: None,
                        youngs_modulus: None,
                        poisson_ratio: None,
                        thermal_expansion: None,
                        reference_temperature: 0.0,
                        num_depvars: 0,
                    });
                }
                Keyword::Expansion => {
                    let kwargs = get_keyword_arguments(line);
                    if let Some(mat) = self.materials.last_mut() {
                        if let Some(zero_str) = kwargs.get("ZERO").and_then(Option::as_deref) {
                            mat.reference_temperature =
                                zero_str.parse().map_err(|e| FerrixError::ParseError {
                                    line: self.line_nr,
                                    message: format!("Invalid ZERO value in *EXPANSION: {e}"),
                                })?;
                        }
                    }
                }
                Keyword::InitialConditions => {
                    let kwargs = get_keyword_arguments(line);
                    self.initial_condition_type = kwargs
                        .get("TYPE")
                        .and_then(Option::as_deref)
                        .map(str::to_uppercase);
                }
                Keyword::SolidSection => {
                    let kwargs = get_keyword_arguments(line);
                    let elset =
                        kwargs
                            .get("ELSET")
                            .and_then(Option::as_deref)
                            .ok_or_else(|| FerrixError::ParseError {
                                line: self.line_nr,
                                message: "Elset not found in *SOLID SECTION".into(),
                            })?;
                    let material_name = kwargs
                        .get("MATERIAL")
                        .and_then(Option::as_deref)
                        .ok_or_else(|| FerrixError::ParseError {
                            line: self.line_nr,
                            message: "Material not found in *SOLID SECTION".into(),
                        })?;

                    let material_index = self
                        .materials
                        .iter()
                        .position(|m| m.name == material_name)
                        .ok_or_else(|| FerrixError::MaterialNotFound(material_name.to_string()))?;

                    if let Some(element_ids) = self.project.mesh.element_sets.get(elset) {
                        for &element_id in element_ids {
                            self.project
                                .element_materials
                                .insert(element_id, material_index);
                        }
                    } else {
                        return Err(FerrixError::ElsetNotFound(elset.to_string()));
                    }
                }
                Keyword::Step => {
                    if self.step_counter == 0 {
                        // Inherit initial state from project
                        self.step_loads.clone_from(&self.project.initial_loads);
                        self.step_bcs.clone_from(&self.project.initial_bcs);
                    }
                    self.step_counter += 1;
                    // If a step was already active, push it
                    if let Some(mut step) = self.current_step.take() {
                        step.loads.clone_from(&self.step_loads);
                        step.bcs.clone_from(&self.step_bcs);
                        self.project.steps.push(Step::StaticStep(step));
                    }

                    let step_args = get_keyword_arguments(line);
                    let nlgeom = step_args.contains_key("NLGEOM");

                    let max_iterations: usize = step_args
                        .get("INC")
                        .and_then(Option::as_deref)
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(10000);

                    if let Some((_next_nr, next_line)) = lines.peek() {
                        if next_line.starts_with("*STATIC") {
                            let step_kwargs = get_keyword_arguments(next_line);
                            let solver = match step_kwargs.get("SOLVER").and_then(Option::as_deref)
                            {
                                Some(solver_str) => match solver_str {
                                    "DIRECT" => SolverType::Direct,
                                    "ITERATIVE" => SolverType::Iterative,
                                    _ => {
                                        return Err(FerrixError::UnknownSolver(
                                            solver_str.to_string(),
                                        ));
                                    }
                                },
                                None => SolverType::Default,
                            };

                            // Consume the *STATIC line
                            lines.next();

                            let mut increment: IncrementData = IncrementData {
                                max_iterations,
                                ..Default::default()
                            };

                            if let Some((_data_nr, data_line)) = lines.peek() {
                                if !data_line.starts_with('*') {
                                    // Increment data was supplied in INP
                                    let args = get_positional_arguments(data_line);
                                    if args.len() >= 4 {
                                        increment.initial_time_increment = args[0]
                                            .parse()
                                            .map_err(|e| FerrixError::ParseError {
                                                line: self.line_nr,
                                                message: format!(
                                                    "Invalid initial time increment: {e}"
                                                ),
                                            })?;
                                        increment.time_period = args[1].parse().map_err(|e| {
                                            FerrixError::ParseError {
                                                line: self.line_nr,
                                                message: format!("Invalid time period: {e}"),
                                            }
                                        })?;
                                        increment.min_time_increment =
                                            args[2].parse().map_err(|e| {
                                                FerrixError::ParseError {
                                                    line: self.line_nr,
                                                    message: format!(
                                                        "Invalid min time increment: {e}"
                                                    ),
                                                }
                                            })?;
                                        increment.max_time_increment =
                                            args[3].parse().map_err(|e| {
                                                FerrixError::ParseError {
                                                    line: self.line_nr,
                                                    message: format!(
                                                        "Invalid max time increment: {e}"
                                                    ),
                                                }
                                            })?;
                                    }
                                    lines.next(); // Consume data line
                                }
                            }

                            self.current_step = Some(StaticStep {
                                solver,
                                increment_data: increment,
                                loads: Vec::new(),
                                bcs: Vec::new(),
                                nlgeom,
                            });
                        }
                    }
                }
                Keyword::EndStep => {
                    if let Some(mut step) = self.current_step.take() {
                        step.loads.clone_from(&self.step_loads);
                        step.bcs.clone_from(&self.step_bcs);
                        self.project.steps.push(Step::StaticStep(step));
                    }
                }
                Keyword::Amplitude => {
                    let kwargs = get_keyword_arguments(line);
                    let name = kwargs
                        .get("NAME")
                        .and_then(Option::as_deref)
                        .ok_or_else(|| FerrixError::ParseError {
                            line: self.line_nr,
                            message: "Amplitude card is missing a Name= argument".into(),
                        })?;
                    let total_time = match kwargs.get("TIME").and_then(Option::as_deref) {
                        Some(val) => val == "TOTAL TIME",
                        None => false,
                    };
                    let shift_x: f64 = kwargs
                        .get("SHIFTX")
                        .and_then(Option::as_deref)
                        .and_then(|val| val.parse().ok())
                        .unwrap_or_default();
                    let shift_y: f64 = kwargs
                        .get("SHIFTY")
                        .and_then(Option::as_deref)
                        .and_then(|val| val.parse().ok())
                        .unwrap_or_default();

                    let mut t = Vec::new();
                    let mut vals = Vec::new();

                    while let Some((_next_nr, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }
                        let (_, line_content) = lines.next().unwrap();
                        let line_t: Vec<f64> = line_content
                            .split(',')
                            .step_by(2)
                            .map(str::trim)
                            .filter(|s| !s.is_empty())
                            .map(|num| {
                                num.parse().map_err(|e| FerrixError::ParseError {
                                    line: self.line_nr,
                                    message: format!("Invalid amplitude time value: {e}"),
                                })
                            })
                            .collect::<Result<Vec<f64>>>()?;

                        let line_vals: Vec<f64> = line_content
                            .split(',')
                            .skip(1)
                            .step_by(2)
                            .map(str::trim)
                            .filter(|s| !s.is_empty())
                            .map(|num| {
                                num.parse().map_err(|e| FerrixError::ParseError {
                                    line: self.line_nr,
                                    message: format!("Invalid amplitude value: {e}"),
                                })
                            })
                            .collect::<Result<Vec<f64>>>()?;

                        t.extend(line_t);
                        vals.extend(line_vals);
                    }

                    let data = Some(TimeSeries(t, vals));

                    self.project.amplitudes.insert(
                        name.into(),
                        Amplitude {
                            total_time,
                            shift_x,
                            shift_y,
                            data,
                        },
                    );
                }
                Keyword::Boundary => {
                    let kwargs = get_keyword_arguments(line);
                    let amplitude_name = kwargs.get("AMPLITUDE").and_then(Option::as_deref);
                    if kwargs
                        .get("OP")
                        .and_then(Option::as_deref)
                        .is_some_and(|op| op == "NEW")
                    {
                        if self.current_step.is_some() {
                            self.step_bcs.clear();
                        } else {
                            self.project.initial_bcs.clear();
                        }
                    }
                    while let Some((_next_nr, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }

                        if let Some((_, line_content)) = lines.next() {
                            self.parse_boundary(line_content, amplitude_name)?;
                        }
                    }
                }
                Keyword::Cload => {
                    let kwargs = get_keyword_arguments(line);
                    let amplitude_name = kwargs.get("AMPLITUDE").and_then(Option::as_deref);
                    if kwargs
                        .get("OP")
                        .and_then(Option::as_deref)
                        .is_some_and(|op| op == "NEW")
                    {
                        if self.current_step.is_some() {
                            self.step_loads.clear();
                        } else {
                            self.project.initial_loads.clear();
                        }
                    }
                    while let Some((_next_nr, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }

                        if let Some((_, line_content)) = lines.next() {
                            self.parse_cload(line_content, amplitude_name)?;
                        }
                    }
                }
                Keyword::NodeFile => {
                    while let Some((_, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }
                        let (_, line_content) = lines.next().unwrap();
                        self.project
                            .nodal_output
                            .extend(line_content.split(',').map(str::trim).map(str::to_string));
                    }
                }
                Keyword::ElFile => {
                    while let Some((_, next_line)) = lines.peek() {
                        if next_line.starts_with('*') {
                            break;
                        }
                        let (_, line_content) = lines.next().unwrap();
                        self.project
                            .element_output
                            .extend(line_content.split(',').map(str::trim).map(str::to_string));
                    }
                }
                _ => {}
            }
        } else {
            return Err(FerrixError::UnsupportedKeyword {
                line: self.line_nr,
                keyword: keyword_str,
            });
        }

        Ok(())
    }

    fn parse_data(&mut self, line: &str) -> Result<()> {
        if let Some(keyword) = self.current_keyword {
            match keyword {
                Keyword::Node => self.parse_node(line)?,
                Keyword::Element => self.parse_element(line)?,
                Keyword::Elastic => self.parse_elastic(line)?,
                Keyword::Density => self.parse_density(line)?,
                Keyword::Expansion => self.parse_expansion(line)?,
                Keyword::InitialConditions => self.parse_initial_conditions(line)?,
                Keyword::Depvar => self.parse_depvar(line)?,
                _ => {}
            }
        }
        Ok(())
    }

    fn parse_node(&mut self, line: &str) -> Result<()> {
        if let Some(node) = Node::parse_line(line) {
            self.project.mesh.nodes.insert(node.id, node);
        } else {
            return Err(FerrixError::ParseError {
                line: self.line_nr,
                message: format!("Invalid node definition: {line}"),
            });
        }
        Ok(())
    }

    fn parse_element(&mut self, line: &str) -> Result<()> {
        if let Some(elem_type) = &self.element_type {
            let elem = Element::parse_line(elem_type, line)?;
            let elem_id = elem.get_id();
            if let Some(elset_name) = &self.elset_name {
                if !elset_name.is_empty() {
                    self.project
                        .mesh
                        .element_sets
                        .entry(elset_name.clone())
                        .or_default()
                        .push(elem_id);
                }
            }
            self.project.mesh.elements.insert(elem_id, elem);
        }
        Ok(())
    }

    fn parse_nset(&mut self, line: &str) {
        if let Some(name) = &self.set_name {
            if self.is_generate {
                let parts: Vec<&str> = line.split(',').map(str::trim).collect();
                if parts.len() >= 2 {
                    let start: usize = parts[0].parse().unwrap_or(0);
                    let end: usize = parts[1].parse().unwrap_or(0);
                    let step: usize = if parts.len() > 2 {
                        parts[2].parse().unwrap_or(1)
                    } else {
                        1
                    };
                    let ids: Vec<NodeId> = (start..=end).step_by(step).map(NodeId).collect();
                    self.project
                        .mesh
                        .node_sets
                        .entry(name.clone())
                        .or_default()
                        .extend(ids);
                }
            } else {
                let ids: Vec<NodeId> = line
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .filter_map(|s| s.parse().ok().map(NodeId))
                    .collect();
                self.project
                    .mesh
                    .node_sets
                    .entry(name.clone())
                    .or_default()
                    .extend(ids);
            }
        }
    }

    fn parse_elset(&mut self, line: &str) {
        if let Some(name) = &self.set_name {
            if self.is_generate {
                let parts: Vec<&str> = line.split(',').map(str::trim).collect();
                if parts.len() >= 2 {
                    let start: usize = parts[0].parse().unwrap_or(0);
                    let end: usize = parts[1].parse().unwrap_or(0);
                    let step: usize = if parts.len() > 2 {
                        parts[2].parse().unwrap_or(1)
                    } else {
                        1
                    };
                    let ids: Vec<ElementId> = (start..=end).step_by(step).map(ElementId).collect();
                    self.project
                        .mesh
                        .element_sets
                        .entry(name.clone())
                        .or_default()
                        .extend(ids);
                }
            } else {
                let parts: Vec<ElementId> = line
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .filter_map(|s| s.parse().ok().map(ElementId))
                    .collect();

                if !line.contains(',') && parts.is_empty() {
                    let other_elset_name = line.trim();
                    let other_ids = self
                        .project
                        .mesh
                        .element_sets
                        .get(other_elset_name)
                        .cloned();
                    if let Some(ids_to_add) = other_ids {
                        self.project
                            .mesh
                            .element_sets
                            .entry(name.clone())
                            .or_default()
                            .extend(ids_to_add);
                    }
                } else {
                    self.project
                        .mesh
                        .element_sets
                        .entry(name.clone())
                        .or_default()
                        .extend(parts);
                }
            }
        }
    }

    fn parse_elastic(&mut self, line: &str) -> Result<()> {
        if let Some(material) = self.materials.last_mut() {
            let parts: Vec<f64> = line
                .split(',')
                .map(|s| {
                    s.trim().parse().map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid elastic constant: {e}"),
                    })
                })
                .collect::<Result<Vec<f64>>>()?;

            if parts.len() >= 2 {
                let e = parts[0];
                let nu = parts[1];
                let temp = if parts.len() >= 3 { parts[2] } else { 0.0 };

                if let Some(lut) = &mut material.youngs_modulus {
                    lut.data.push((temp, e));
                    lut.data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                } else {
                    material.youngs_modulus = Some(TemperatureDependentLUT::new(vec![(temp, e)]));
                }

                if let Some(lut) = &mut material.poisson_ratio {
                    lut.data.push((temp, nu));
                    lut.data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                } else {
                    material.poisson_ratio = Some(TemperatureDependentLUT::new(vec![(temp, nu)]));
                }
            }
        }
        Ok(())
    }

    fn parse_density(&mut self, line: &str) -> Result<()> {
        if let Some(material) = self.materials.last_mut() {
            let parts: Vec<f64> = line
                .split(',')
                .map(|s| {
                    s.trim().parse().map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid density value: {e}"),
                    })
                })
                .collect::<Result<Vec<f64>>>()?;

            if !parts.is_empty() {
                let rho = parts[0];
                let temp = if parts.len() >= 2 { parts[1] } else { 0.0 };

                if let Some(lut) = &mut material.density {
                    lut.data.push((temp, rho));
                    lut.data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                } else {
                    material.density = Some(TemperatureDependentLUT::new(vec![(temp, rho)]));
                }
            }
        }
        Ok(())
    }

    fn parse_expansion(&mut self, line: &str) -> Result<()> {
        if let Some(material) = self.materials.last_mut() {
            let parts: Vec<f64> = line
                .split(',')
                .map(|s| {
                    s.trim().parse().map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid expansion constant: {e}"),
                    })
                })
                .collect::<Result<Vec<f64>>>()?;

            if !parts.is_empty() {
                let alpha = parts[0];
                let temp = if parts.len() >= 2 { parts[1] } else { 0.0 };

                if let Some(lut) = &mut material.thermal_expansion {
                    lut.data.push((temp, alpha));
                    lut.data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                } else {
                    material.thermal_expansion =
                        Some(TemperatureDependentLUT::new(vec![(temp, alpha)]));
                }
            }
        }
        Ok(())
    }

    fn parse_initial_conditions(&mut self, line: &str) -> Result<()> {
        if let Some(ic_type) = &self.initial_condition_type {
            if ic_type == "TEMPERATURE" {
                let parts: Vec<&str> = line.split(',').map(str::trim).collect();
                if parts.len() >= 2 {
                    let target = parts[0];
                    let temp: f64 = parts[1].parse().map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid temperature value: {e}"),
                    })?;

                    let resolve_target = |target: &str, mesh: &Mesh| -> Vec<NodeId> {
                        let t = target.trim();
                        if let Some(ids) = mesh.node_sets.get(t) {
                            return ids.clone();
                        }
                        if let Ok(id) = t.parse::<usize>() {
                            return vec![NodeId(id)];
                        }
                        Vec::new()
                    };

                    let target_nodes = resolve_target(target, &self.project.mesh);
                    for node_id in target_nodes {
                        self.project.initial_temperatures.insert(node_id, temp);
                    }
                }
            }
        }
        Ok(())
    }

    fn parse_depvar(&mut self, line: &str) -> Result<()> {
        if let Some(material) = self.materials.last_mut() {
            let parts: Vec<&str> = line.split(',').map(str::trim).collect();
            if !parts.is_empty() {
                let n_sdv: usize = parts[0].parse().map_err(|e| FerrixError::ParseError {
                    line: self.line_nr,
                    message: format!("Invalid DEPVAR value: {e}"),
                })?;
                material.num_depvars = n_sdv;
            }
        }
        Ok(())
    }

    fn parse_cload(&mut self, line: &str, amplitude_name: Option<&str>) -> Result<()> {
        let parts: Vec<&str> = line.split(',').map(str::trim).collect();
        if parts.len() >= 3 {
            let target = parts[0];
            let dof_in: usize = parts[1]
                .trim()
                .parse()
                .map_err(|e| FerrixError::ParseError {
                    line: self.line_nr,
                    message: format!("Invalid DOF: {e}"),
                })?;
            let val: f64 = parts[2]
                .trim()
                .parse()
                .map_err(|e| FerrixError::ParseError {
                    line: self.line_nr,
                    message: format!("Invalid load value: {e}"),
                })?;

            let resolve_target = |target: &str, mesh: &Mesh| -> Vec<NodeId> {
                let t = target.trim();
                if let Some(ids) = mesh.node_sets.get(t) {
                    return ids.clone();
                }
                if let Ok(id) = t.parse::<usize>() {
                    return vec![NodeId(id)];
                }
                Vec::new()
            };

            let target_nodes = resolve_target(target, &self.project.mesh);

            if (1..=3).contains(&dof_in) {
                for node_id in target_nodes {
                    let load = Load::new(
                        node_id,
                        dof_in - 1,
                        val,
                        match amplitude_name {
                            Some(name) => self.project.amplitudes.get(name).cloned(),
                            None => None,
                        },
                        if self.current_step.is_some() {
                            self.step_counter
                        } else {
                            0
                        },
                    );
                    if self.current_step.is_some() {
                        self.step_loads.push(load);
                    } else {
                        self.project.initial_loads.push(load);
                    }
                    self.load_id_counter += 1;
                }
            }
        }
        Ok(())
    }

    fn parse_boundary(&mut self, line: &str, amplitude_name: Option<&str>) -> Result<()> {
        let parts: Vec<&str> = line.split(',').map(str::trim).collect();
        if parts.len() >= 2 {
            let target = parts[0];
            let first_dof: usize =
                parts[1]
                    .trim()
                    .parse()
                    .map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid first DOF: {e}"),
                    })?;
            let last_dof: usize = if parts.len() > 2 && !parts[2].trim().is_empty() {
                parts[2]
                    .trim()
                    .parse()
                    .map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid last DOF: {e}"),
                    })?
            } else {
                first_dof
            };
            let val: f64 = if parts.len() > 3 {
                parts[3]
                    .trim()
                    .parse()
                    .map_err(|e| FerrixError::ParseError {
                        line: self.line_nr,
                        message: format!("Invalid boundary value: {e}"),
                    })?
            } else {
                0.0
            };

            let resolve_target = |target: &str, mesh: &Mesh| -> Vec<NodeId> {
                let t = target.trim();
                if let Some(ids) = mesh.node_sets.get(t) {
                    return ids.clone();
                }
                if let Ok(id) = t.parse::<usize>() {
                    return vec![NodeId(id)];
                }
                Vec::new()
            };

            let target_nodes = resolve_target(target, &self.project.mesh);

            for node_id in target_nodes {
                for dof_in in first_dof..=last_dof {
                    if (1..=3).contains(&dof_in) {
                        let bc = BoundaryCondition::new(
                            node_id,
                            dof_in - 1,
                            val,
                            match amplitude_name {
                                Some(name) => self.project.amplitudes.get(name).cloned(),
                                None => None,
                            },
                            if self.current_step.is_some() {
                                self.step_counter
                            } else {
                                0
                            },
                        );
                        if self.current_step.is_some() {
                            self.step_bcs.push(bc);
                        } else {
                            self.project.initial_bcs.push(bc);
                        }
                        self.bc_id_counter += 1;
                    }
                }
            }
        }
        Ok(())
    }
}

fn get_keyword_arguments(line: &str) -> HashMap<&str, Option<&str>> {
    line.split(',')
        .skip(1) // Skip the keyword itself (*STEP, etc)
        .map(|s| {
            if let Some((k, v)) = s.split_once('=') {
                (k.trim(), Some(v.trim()))
            } else {
                (s.trim(), None)
            }
        })
        .collect()
}

fn get_positional_arguments(line: &str) -> Vec<&str> {
    line.split(',').map(str::trim).collect()
}

/// This preprocess includes:
/// - removing leading and trailing whitespaces
/// - removing comments
/// - making all text uppercase
/// - merging lines that belong together (`,` at the end of line)
/// - removes empty lines
#[must_use]
pub fn preprocess_inp(input_file: &str) -> String {
    let mut preprocessed = String::new();
    let mut current_line = String::new();

    for line in input_file.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("**") {
            continue;
        }

        let upper_line = line.to_uppercase();
        if !current_line.is_empty() {
            current_line.push(' ');
        }
        current_line.push_str(&upper_line);

        if !upper_line.ends_with(',') {
            preprocessed.push_str(&current_line);
            preprocessed.push('\n');
            current_line.clear();
        }
    }

    if !current_line.is_empty() {
        preprocessed.push_str(&current_line);
        preprocessed.push('\n');
    }

    substitute_ccx_solvers(&preprocessed)
}

use aho_corasick::AhoCorasick;

/// Maps the common `CalculiX` Solvers to supported `Ferrix` solvers.
///
/// # Panics
/// Panics if the internal `AhoCorasick` engine fails to initialize.
#[must_use]
pub fn substitute_ccx_solvers(input_file: &str) -> String {
    let patterns = &[
        "ITERATIVE CHOLESKY",
        "ITERATIVE SCALING",
        "PASTIX",
        "PARDISO",
        "SPOOLES",
    ];
    let replaces = &["ITERATIVE", "ITERATIVE", "DIRECT", "DIRECT", "DIRECT"];

    let ac = AhoCorasick::new(patterns).expect("Failed to initialize AhoCorasick");
    ac.replace_all(input_file, replaces)
}

#[cfg(test)]
mod tests {
    use itertools::assert_equal;

    use super::*;

    #[test]
    fn test_preprocess_inp() {
        let inp = "**comment\n word  \n \t*keyword\n123.4\n4, 5, 6,\n7, 8, 9";
        assert_eq!(
            preprocess_inp(inp),
            "WORD\n*KEYWORD\n123.4\n4, 5, 6, 7, 8, 9\n"
        );
    }

    #[test]
    fn keyword_args() {
        let line = "*STEP , INC =  100";
        let args = get_keyword_arguments(line);
        assert!(args["INC"] == Some("100"));
    }

    #[test]
    fn positional_args() {
        let line = "0.1, 1, 1E-05, 0.2";
        let args = get_positional_arguments(line);
        assert_equal(args, vec!["0.1", "1", "1E-05", "0.2"]);
    }

    #[test]
    fn test_multiple_boundary_lines() {
        let inp = "
*NODE
1, 0, 0, 0
2, 1, 0, 0
*BOUNDARY
1, 1, 1, 0.1
1, 2, 2, 0.2
*STEP
*STATIC
*BOUNDARY
2, 1, 1, 0.3
2, 2, 2, 0.4
*END STEP
";
        let inp_file = InpFile::new(inp);
        let parser = Parser::new(&inp_file);
        let project = parser.parse().unwrap();

        assert_eq!(project.initial_bcs.len(), 2, "Should have 2 initial BCs");
        assert_eq!(project.steps.len(), 1);
        let Step::StaticStep(step) = &project.steps[0];
        assert_eq!(
            step.bcs.len(),
            4,
            "Should have 4 step BCs (2 inherited + 2 new)"
        );
    }

    #[test]
    fn test_multiple_amplitude_lines() {
        let inp = "
        *AMPLITUDE, NAME=TABULAR-1
        0, 0, 1, 100,
        2, 200, 3, 300
        ";
        let inp_file = InpFile::new(inp);
        let parser = Parser::new(&inp_file);
        let project = parser.parse().unwrap();

        let amp = project.amplitudes.get("TABULAR-1").unwrap();
        if let Some(TimeSeries(t, vals)) = &amp.data {
            assert_eq!(t.len(), 4);
            assert_eq!(vals.len(), 4);
            assert_eq!(t, &vec![0.0, 1.0, 2.0, 3.0]);
            assert_eq!(vals, &vec![0.0, 100.0, 200.0, 300.0]);
        } else {
            panic!("Amplitude data should be present");
        }
    }
}
