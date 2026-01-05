use std::{
    collections::{HashMap, HashSet},
    error::Error,
};

use nalgebra::{Point3, center};

use crate::solver::{
    mesh_lib::{
        elements::element::{Element, ElementType, Face},
        node::Node,
    },
    project::{InpParsingError, InpSection},
};

/// Contains all Node and Element Data
#[derive(Debug)]
pub struct Mesh {
    pub nodes: HashMap<usize, Node>,
    pub elements: Vec<Element>,
    pub wireframe_lines: Vec<(Point3<f64>, Point3<f64>)>,
}

impl Mesh {
    #[allow(clippy::match_wildcard_for_single_variants)]
    pub fn from_sections(
        input_file: &str,
        sections: &Vec<InpSection>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut elements = Vec::new();
        let mut nodes: Option<Vec<Node>> = None;
        for sec in sections {
            match sec {
                InpSection::Node(nr) => {
                    nodes = Some(
                        input_file
                            .lines()
                            .skip(*nr + 1)
                            .map_while(Node::parse_line)
                            .collect(),
                    );
                }
                InpSection::Element(nr) => {
                    let elem_type = Element::parse_type_str_from_line(
                        input_file
                            .lines()
                            .nth(*nr)
                            .expect("The line number is outside the file, aborting"),
                    )?;
                    elements.extend(
                        input_file
                            .lines()
                            .skip(nr + 1)
                            .take_while(|line| {
                                line.chars()
                                    .nth(0)
                                    .expect("There are no empty lines after preprocessing")
                                    .is_numeric()
                            })
                            .map(|line| Element::parse_line(&elem_type, line)),
                    );
                }

                _ => {}
            }
        }
        let mut node_hash: HashMap<usize, Node> = HashMap::new();
        for node in
            nodes.ok_or("The input file does not contain a *NODE card. Analysis is aborted.")?
        {
            node_hash.insert(node.id, node);
        }

        let mut mesh = Self {
            nodes: node_hash,
            elements,
            wireframe_lines: Vec::new(),
        };
        mesh.precompute_wireframe();
        Ok(mesh)
    }

    /// Counts all elements by their respective type.
    pub fn count_by_type(&self) -> HashMap<ElementType, u32> {
        let mut elem_count: HashMap<ElementType, u32> = HashMap::new();
        for elem in &self.elements {
            let elem_type: ElementType = elem.into();
            *elem_count.entry(elem_type).or_insert(0) += 1;
        }
        elem_count
    }

    /// Compute the median Position of all points
    pub fn get_center(&self) -> Point3<f64> {
        if self.nodes.is_empty() {
            return Point3::origin();
        }

        let init_min = Point3::new(f64::MAX, f64::MAX, f64::MAX);
        let init_max = Point3::new(f64::MIN, f64::MIN, f64::MIN);

        let (min, max) = self
            .nodes
            .iter()
            .fold((init_min, init_max), |(min, max), (_, node)| {
                (
                    Point3::new(min.x.min(node.x), min.y.min(node.y), min.z.min(node.z)),
                    Point3::new(max.x.max(node.x), max.y.max(node.y), max.z.max(node.z)),
                )
            });

        center(&min, &max)
    }

    pub fn precompute_wireframe(&mut self) {
        let mut face_counts: HashMap<Face, u8> = HashMap::new();

        // 1. Zähle Vorkommen aller Flächen
        for elem in &self.elements {
            for face in elem.get_faces() {
                // Wir nutzen saturating_add, da uns alles >= 2 egal ist (es ist internal)
                face_counts
                    .entry(face)
                    .and_modify(|c| *c = c.saturating_add(1))
                    .or_insert(1);
            }
        }

        // 2. Filtere nur Flächen mit Count == 1 (Außenhaut)
        // 3. Sammle Kanten dieser Flächen in ein HashSet (zur Deduplizierung)
        let mut unique_edges: HashSet<(usize, usize)> = HashSet::new();

        for (face, count) in face_counts {
            if count == 1 {
                let ids = match face {
                    Face::Tri(n) => vec![(n[0], n[1]), (n[1], n[2]), (n[2], n[0])],
                    Face::Quad(n) => vec![(n[0], n[1]), (n[1], n[2]), (n[2], n[3]), (n[3], n[0])],
                };

                for (a, b) in ids {
                    // Sortiere Kanten-IDs (min, max), damit (1,2) == (2,1)
                    let edge = if a < b { (a, b) } else { (b, a) };
                    unique_edges.insert(edge);
                }
            }
        }

        // 4. Löse IDs in Koordinaten auf und speichere im Cache
        self.wireframe_lines = unique_edges
            .iter()
            .filter_map(|(id_a, id_b)| {
                let n_a = self.nodes.get(id_a)?;
                let n_b = self.nodes.get(id_b)?;
                Some((
                    Point3::new(n_a.x, n_a.y, n_a.z),
                    Point3::new(n_b.x, n_b.y, n_b.z),
                ))
            })
            .collect();
    }
}
