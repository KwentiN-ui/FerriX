use crate::solver::ids::{BoundaryConditionId, ElementId, LoadId, NodeId};
use crate::solver::inp::InpFile;
use crate::solver::material::Material;
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
use strum_macros::EnumString;

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
}

pub struct Parser<'a> {
    input: &'a InpFile,
    project: Project,
    current_keyword: Option<Keyword>,
    line_nr: usize,
    // Parser state
    element_type: Option<String>,
    elset_name: Option<String>,
    set_name: Option<String>,
    is_generate: bool,
    load_id_counter: usize,
    bc_id_counter: usize,
}

impl<'a> Parser<'a> {
    pub fn new(input: &'a InpFile) -> Self {
        Self {
            input,
            project: Project::new(),
            current_keyword: None,
            line_nr: 0,
            element_type: None,
            elset_name: None,
            set_name: None,
            is_generate: false,
            load_id_counter: 0,
            bc_id_counter: 0,
        }
    }

    pub fn parse(mut self) -> Result<Project, String> {
        let mut lines = self.input.lines().enumerate().peekable();
        while let Some((line_nr, line_content)) = lines.next() {
            self.line_nr = line_nr;

            if line_content.starts_with('*') {
                self.parse_keyword(line_content, &mut lines)?;
                continue;
            }

            let line_content = line_content.trim();
            if line_content.is_empty() {
                continue;
            }

            if self.current_keyword.is_some() {
                self.parse_data(line_content);
            }
        }

        // Post-parsing steps, e.g. building node mappings
        self.project.mesh.build_node_mappings();

        Ok(self.project)
    }

    fn parse_keyword(
        &mut self,
        line: &str,
        lines: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines>>,
    ) -> Result<(), String> {
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
                    self.element_type = Element::parse_type_str_from_line(line).ok();
                    self.elset_name = parts
                        .iter()
                        .find(|p| p.starts_with("ELSET="))
                        .map(|p| p.split('=').nth(1).unwrap_or("").trim().to_string());
                }
                Keyword::Nset => {
                    self.set_name = parts
                        .iter()
                        .find(|p| p.starts_with("NSET="))
                        .map(|p| p.split('=').nth(1).unwrap_or("").trim().to_string());
                    self.is_generate = parts.iter().any(|p| p.to_uppercase() == "GENERATE");
                }
                Keyword::Elset => {
                    self.set_name = parts
                        .iter()
                        .find(|p| p.starts_with("ELSET="))
                        .map(|p| p.split('=').nth(1).unwrap_or("").trim().to_string());
                    self.is_generate = parts.iter().any(|p| p.to_uppercase() == "GENERATE");
                }
                Keyword::Material => {
                    let name = parts
                        .iter()
                        .find(|p| p.starts_with("NAME="))
                        .map(|p| p.split('=').nth(1).unwrap_or("").trim().to_string())
                        .ok_or("Material name not found")?;
                    self.project.materials.push(Material {
                        name,
                        density: None,
                        elastic: None,
                    });
                }
                Keyword::SolidSection => {
                    let elset = parts
                        .iter()
                        .find(|p| p.starts_with("ELSET="))
                        .map(|p| p.split('=').nth(1).unwrap_or("").trim().to_string())
                        .ok_or("Elset not found in *SOLID SECTION")?;
                    let material_name = parts
                        .iter()
                        .find(|p| p.starts_with("MATERIAL="))
                        .map(|p| p.split('=').nth(1).unwrap_or("").trim().to_string())
                        .ok_or("Material not found in *SOLID SECTION")?;

                    let material_index = self
                        .project
                        .materials
                        .iter()
                        .position(|m| m.name == material_name)
                        .ok_or(format!(
                            "Material {material_name} not found for *SOLID SECTION"
                        ))?;

                    if let Some(element_ids) = self.project.mesh.element_sets.get(&elset) {
                        for &element_id in element_ids {
                            self.project
                                .element_materials
                                .insert(element_id, material_index);
                        }
                    } else {
                        return Err(format!("Elset {elset} not found for *SOLID SECTION"));
                    }
                }
                Keyword::Step => {
                    let max_increments: usize = parse_keyword_arguments(line)
                        .get("INC")
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(100);

                    if let Some((_next_nr, next_line)) = lines.peek() {
                        if next_line.starts_with("*STATIC") {
                            let step_kwargs = parse_keyword_arguments(next_line);
                            let solver = match step_kwargs.get("SOLVER") {
                                Some(solver_str) => match *solver_str {
                                    "DIRECT" => SolverType::Direct,
                                    "ITERATIVE" => SolverType::Iterative,
                                    _ => panic!("Unknown solver type: {}", solver_str),
                                },
                                None => SolverType::Default,
                            };

                            // Consume the *STATIC line
                            lines.next();

                            if let Some((_data_nr, data_line)) = lines.next() {
                                let time_data: Vec<f64> = data_line
                                    .split(',')
                                    .filter_map(|s| s.trim().parse().ok())
                                    .collect();

                                if time_data.len() >= 4 {
                                    let static_step = StaticStep {
                                        solver,
                                        max_increments,
                                        initial_time_increment: time_data[0],
                                        time_period: time_data[1],
                                        min_time_increment: time_data[2],
                                        max_time_increment: time_data[3],
                                    };
                                    self.project.steps.push(Step::StaticStep(static_step));
                                } else {
                                    // Handle error: not enough time data
                                    return Err(
                                        "Not enough values for *STATIC time data".to_string()
                                    );
                                }
                            } else {
                                // Handle error: no data line after *STATIC
                                return Err("Expected data line after *STATIC".to_string());
                            }
                        }
                    }
                }
                _ => {}
            }
        } else {
            self.current_keyword = None;
        }

        Ok(())
    }

    fn parse_data(&mut self, line: &str) {
        if let Some(keyword) = self.current_keyword {
            match keyword {
                Keyword::Node => self.parse_node(line),
                Keyword::Element => self.parse_element(line),
                Keyword::Nset => self.parse_nset(line),
                Keyword::Elset => self.parse_elset(line),
                Keyword::Elastic => self.parse_elastic(line),
                Keyword::Cload => self.parse_cload(line),
                Keyword::Boundary => self.parse_boundary(line),
                Keyword::NodeFile => {
                    self.project
                        .nodal_output
                        .extend(line.split(',').map(str::trim).map(str::to_string));
                }
                Keyword::ElFile => {
                    self.project
                        .element_output
                        .extend(line.split(',').map(str::trim).map(str::to_string));
                }
                _ => {} // For now
            }
        }
    }

    fn parse_node(&mut self, line: &str) {
        if let Some(node) = Node::parse_line(line) {
            self.project.mesh.nodes.insert(node.id, node);
        }
    }

    fn parse_element(&mut self, line: &str) {
        if let Some(elem_type) = &self.element_type {
            let elem = Element::parse_line(elem_type, line);
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
    }

    fn parse_nset(&mut self, line: &str) {
        if let Some(name) = &self.set_name {
            if self.is_generate {
                // TODO
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
                // TODO
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

    fn parse_elastic(&mut self, line: &str) {
        if let Some(material) = self.project.materials.last_mut() {
            let parts: Vec<f64> = line
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if parts.len() >= 2 {
                material.elastic = Some((parts[0], parts[1]));
            }
        }
    }

    fn parse_cload(&mut self, line: &str) {
        let parts: Vec<&str> = line.split(',').map(str::trim).collect();
        if parts.len() >= 3 {
            let target = parts[0];
            let dof_in: usize = parts[1].trim().parse().unwrap_or(0);
            let val: f64 = parts[2].trim().parse().unwrap_or(0.0);

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
                    self.project.loads.push(Load {
                        id: LoadId(self.load_id_counter),
                        node_id,
                        dof: dof_in - 1,
                        value: val,
                    });
                    self.load_id_counter += 1;
                }
            }
        }
    }

    fn parse_boundary(&mut self, line: &str) {
        let parts: Vec<&str> = line.split(',').map(str::trim).collect();
        if parts.len() >= 2 {
            let target = parts[0];
            let first_dof: usize = parts[1].trim().parse().unwrap_or(0);
            let last_dof: usize = if parts.len() > 2 && !parts[2].trim().is_empty() {
                parts[2].trim().parse().unwrap_or(first_dof)
            } else {
                first_dof
            };
            let val: f64 = if parts.len() > 3 {
                parts[3].trim().parse().unwrap_or(0.0)
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
                        self.project.bcs.push(BoundaryCondition {
                            id: BoundaryConditionId(self.bc_id_counter),
                            node_id,
                            dof: dof_in - 1,
                            value: val,
                        });
                        self.bc_id_counter += 1;
                    }
                }
            }
        }
    }
}

fn parse_keyword_arguments(line: &str) -> HashMap<&str, &str> {
    line.split(',')
        .filter_map(|s| s.split_once('='))
        .map(|(k, v)| (k.trim().into(), v.trim().into()))
        .collect()
}

fn parse_postional_arguments(line: &str) -> Vec<&str> {
    line.split(',').map(str::trim).collect()
}

/// This preprocess includes:
/// - removing leading and trailing whitespaces
/// - removing comments
/// - making all text uppercase
/// - merging lines that belong together (`,` at the end of line)
/// - removes empty lines
pub fn preprocess_inp(input_file: &str) -> String {
    let preprocessed: String = input_file
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(|line| line.chars().map(|c| c.to_uppercase().to_string()).collect())
        // merge lines that end with `,`
        .map(|line: String| {
            if line.ends_with(',') {
                line + " "
            } else {
                line + "\n"
            }
        })
        // remove comments
        .filter(|line| !line.starts_with("**"))
        .collect();
    substitute_ccx_solvers(&preprocessed)
}

use aho_corasick::AhoCorasick;

/// Maps the common `CalculiX` Solvers to supported `Ferrix` solvers.
pub fn substitute_ccx_solvers(input_file: &str) -> String {
    let patterns = &[
        "ITERATIVE CHOLESKY",
        "ITERATIVE SCALING",
        "PASTIX",
        "PARDISO",
        "SPOOLES",
    ];
    let replaces = &["ITERATIVE", "ITERATIVE", "DIRECT", "DIRECT", "DIRECT"];

    let ac = AhoCorasick::new(patterns).unwrap();
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
        let args = parse_keyword_arguments(line);
        assert!(args["INC"] == "100".to_string());
    }

    #[test]
    fn positional_args() {
        let line = "0.1, 1, 1E-05, 0.2";
        let args = parse_postional_arguments(line);
        assert_equal(
            args,
            vec![
                "0.1".to_string(),
                "1".to_string(),
                "1E-05".to_string(),
                "0.2".to_string(),
            ],
        );
    }
}
