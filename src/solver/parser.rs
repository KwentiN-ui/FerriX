use crate::solver::inp::InpFile;
use crate::solver::material::Material;
use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::mesh_lib::node::Node;
use crate::solver::project::Project;
use crate::solver::step::steps::Step;
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
        }
    }

    pub fn parse(mut self) -> Result<Project, String> {
        let mut lines = self.input.0.lines().enumerate().peekable();
        while let Some((line_nr, line_content)) = lines.next() {
            self.line_nr = line_nr;
            if line_content.starts_with("**") {
                continue;
            }

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
        
        // Post-parsing steps, e.g. building node mappings
        self.project.mesh.build_node_mappings();

        Ok(self.project)
    }

    fn parse_keyword(&mut self, line: &str, lines: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines>>) -> Result<(), String> {
        let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        let keyword_str = parts[0].strip_prefix('*').unwrap_or("").trim().to_uppercase();

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
                        .ok_or(format!("Material {} not found for *SOLID SECTION", material_name))?;
                    
                    if let Some(element_ids) = self.project.mesh.element_sets.get(&elset) {
                        for &element_id in element_ids {
                            self.project.element_materials.insert(element_id, material_index);
                        }
                    } else {
                        return Err(format!("Elset {} not found for *SOLID SECTION", elset));
                    }
                }
                Keyword::Step => {
                    if let Some((next_nr, next_line)) = lines.peek() {
                        if next_line.starts_with("*STATIC") {
                            self.project.steps.push(Step::StaticStep(*next_nr));
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

    fn parse_data(&mut self, line: &str) -> Result<(), String> {
        if let Some(keyword) = self.current_keyword {
            match keyword {
                Keyword::Node => self.parse_node(line)?,
                Keyword::Element => self.parse_element(line)?,
                Keyword::Nset => self.parse_nset(line)?,
                Keyword::Elset => self.parse_elset(line)?,
                Keyword::Elastic => self.parse_elastic(line)?,
                _ => {} // For now
            }
        }
        Ok(())
    }

    fn parse_node(&mut self, line: &str) -> Result<(), String> {
        if let Some(node) = Node::parse_line(line) {
            self.project.mesh.nodes.insert(node.id, node);
        }
        Ok(())
    }

    fn parse_element(&mut self, line: &str) -> Result<(), String> {
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
        Ok(())
    }

    fn parse_nset(&mut self, line: &str) -> Result<(), String> {
        if let Some(name) = &self.set_name {
            if self.is_generate {
                // TODO
            }
            else {
                let ids: Vec<usize> = line
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .filter_map(|s| s.parse().ok())
                    .collect();
                self.project.mesh.node_sets.entry(name.clone()).or_default().extend(ids);
            }
        }
        Ok(())
    }

    fn parse_elset(&mut self, line: &str) -> Result<(), String> {
        if let Some(name) = &self.set_name {
            if self.is_generate {
                // TODO
            }
            else {
                let parts: Vec<usize> = line
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .filter_map(|s| s.parse::<usize>().ok())
                    .collect();

                if !line.contains(',') && parts.is_empty() {
                    let other_elset_name = line.trim();
                    let other_ids = self.project.mesh.element_sets.get(other_elset_name).cloned();
                    if let Some(ids_to_add) = other_ids {
                        self.project.mesh.element_sets.entry(name.clone()).or_default().extend(ids_to_add);
                    }
                } else {
                    self.project.mesh.element_sets.entry(name.clone()).or_default().extend(parts);
                }
            }
        }
        Ok(())
    }

    fn parse_elastic(&mut self, line: &str) -> Result<(), String> {
        if let Some(material) = self.project.materials.last_mut() {
            let parts: Vec<f64> = line
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if parts.len() >= 2 {
                material.elastic = Some((parts[0], parts[1]));
            }
        }
        Ok(())
    }
}