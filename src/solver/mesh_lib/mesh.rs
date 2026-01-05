use std::{
    collections::{HashMap, HashSet},
    error::Error,
};

use nalgebra::{Point3, Unit, Vector3, center};

use crate::solver::{
    inp::InpFile,
    mesh_lib::{
        elements::element::{Element, ElementType, Face},
        node::Node,
    },
    project::{InpParsingError, InpSection},
};

/// Contains all Node and Element Data
#[derive(Debug, Clone)]
pub struct Mesh {
    pub nodes: HashMap<usize, Node>,
    pub elements: HashMap<usize, Element>,
    pub wireframe_lines: Vec<(Point3<f64>, Point3<f64>)>,
    pub node_sets: HashMap<String, Vec<usize>>,

    /// Map: Node-ID (from INP) -> Matrix-Index (0..N)
    pub node_id_to_index: HashMap<usize, usize>,
    /// Map: Matrix-Index (0..N) -> Node-ID (from INP)
    pub index_to_node_id: Vec<usize>,
}

impl Mesh {
    pub fn new() -> Self {
        Self {
            nodes: HashMap::new(),
            elements: HashMap::new(),
            wireframe_lines: Vec::new(),
            node_id_to_index: HashMap::new(),
            index_to_node_id: Vec::new(),
            node_sets: HashMap::new(),
        }
    }
    pub fn from_sections(
        input_file: &InpFile,
        sections: &Vec<InpSection>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut elements = Vec::new();
        let mut nodes: Option<Vec<Node>> = None;
        let mut node_sets: HashMap<String, Vec<usize>> = HashMap::new();

        for sec in sections {
            match sec {
                InpSection::Node(nr) => {
                    nodes = Some(
                        input_file
                            .0
                            .lines()
                            .skip(*nr + 1)
                            .map_while(Node::parse_line)
                            .collect(),
                    );
                }
                InpSection::Element(nr) => {
                    let elem_type = Element::parse_type_str_from_line(
                        input_file
                            .0
                            .lines()
                            .nth(*nr)
                            .expect("Line number outside file"),
                    )?;
                    elements.extend(
                        input_file
                            .0
                            .lines()
                            .skip(nr + 1)
                            .take_while(|line| {
                                line.chars().nth(0).expect("No empty lines").is_numeric()
                            })
                            .map(|line| Element::parse_line(&elem_type, line)),
                    );
                }
                // Parse Node Sets
                // Assumes InpSection::Nset(line_index, name, is_generate)
                InpSection::Nset(nr, name, is_generate) => {
                    let mut ids = Vec::new();

                    // Read lines until next keyword (*)
                    let data_lines = input_file
                        .0
                        .lines()
                        .skip(nr + 1)
                        .take_while(|line| !line.trim().starts_with('*'));

                    if *is_generate {
                        // Format: Start, End, Increment
                        for line in data_lines {
                            let nums: Vec<usize> = line
                                .split(',')
                                .map(|s| s.trim())
                                .filter(|s| !s.is_empty())
                                .filter_map(|s| s.parse().ok())
                                .collect();

                            if nums.len() >= 2 {
                                // At least Start, End
                                let start = nums[0];
                                let end = nums[1];
                                let step = *nums.get(2).unwrap_or(&1); // Default step 1

                                for id in (start..=end).step_by(step) {
                                    ids.push(id);
                                }
                            }
                        }
                    } else {
                        // Explicit List: 1, 2, 3, 4...
                        for line in data_lines {
                            let line_ids = line
                                .split(',')
                                .map(|s| s.trim())
                                .filter(|s| !s.is_empty())
                                .filter_map(|s| s.parse::<usize>().ok());
                            ids.extend(line_ids);
                        }
                    }
                    node_sets.insert(name.clone(), ids);
                }

                _ => {}
            }
        }

        let mut node_hash: HashMap<usize, Node> = HashMap::new();
        for node in nodes.ok_or("Input file missing *NODE card.")? {
            node_hash.insert(node.id, node);
        }

        let mut elem_hash: HashMap<usize, Element> = HashMap::new();
        for elem in elements {
            elem_hash.insert(elem.get_id(), elem);
        }

        // 3. Construct Mesh with populated sets
        let mut mesh = Self {
            nodes: node_hash,
            elements: elem_hash,
            node_sets, // Pass the filled map
            wireframe_lines: Vec::new(),
            node_id_to_index: HashMap::new(),
            index_to_node_id: Vec::new(),
        };

        mesh.precompute_wireframe();
        mesh.build_node_mappings();

        Ok(mesh)
    }

    pub fn build_node_mappings(&mut self) {
        let mut sorted_ids: Vec<usize> = self.nodes.keys().copied().collect();
        sorted_ids.sort_unstable();

        // Reverse Mapping (Index -> ID)
        self.index_to_node_id.clone_from(&sorted_ids);

        // Forward Mapping (ID -> Index)
        self.node_id_to_index = sorted_ids
            .into_iter()
            .enumerate()
            .map(|(idx, id)| (id, idx))
            .collect();
    }

    /// Fetches the matrix-index for a given Node-ID (from INP file)
    pub fn get_index_for_node_id(&self, id: usize) -> Option<usize> {
        self.node_id_to_index.get(&id).copied()
    }

    /// Resolves a string (from *BOUNDARY) to a list of Node IDs.
    /// The string can be a direct number ("105") or a Set Name ("FIXPOINTS").
    pub fn resolve_node_ids(&self, identifier: &str) -> Vec<usize> {
        let identifier = identifier.trim();

        // 1. Try to find a Node Set with this name
        if let Some(set_ids) = self.node_sets.get(identifier) {
            return set_ids.clone();
        }

        // 2. Try to parse as single Node ID
        if let Ok(id) = identifier.parse::<usize>() {
            // Optional: Check if node exists? For now just return it.
            return vec![id];
        }

        // 3. Not found
        Vec::new()
    }

    /// Counts all elements by their respective type.
    pub fn count_by_type(&self) -> HashMap<ElementType, u32> {
        let mut elem_count: HashMap<ElementType, u32> = HashMap::new();
        for elem in self.elements.values() {
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

    /// This method creates a wireframe for hard edges of the model for the tui mesh preview window.
    pub fn precompute_wireframe(&mut self) {
        // Surface Extraction
        let mut face_counts: HashMap<Face, u8> = HashMap::new();
        for elem in self.elements.values() {
            for face in elem.get_faces() {
                face_counts
                    .entry(face)
                    .and_modify(|c| *c = c.saturating_add(1))
                    .or_insert(1);
            }
        }

        // Map: Edge -> List of normals
        // Edge Key is (min_id, max_id)
        let mut edge_normals: HashMap<(usize, usize), Vec<Vector3<f64>>> = HashMap::new();

        // iterate over external faces (Count == 1)
        for (face, count) in face_counts {
            if count == 1 {
                if let Some(normal) = self.compute_face_normal(&face) {
                    // get all edges of the face
                    let ids = match face {
                        Face::Tri(n) => vec![(n[0], n[1]), (n[1], n[2]), (n[2], n[0])],
                        Face::Quad(n) => {
                            vec![(n[0], n[1]), (n[1], n[2]), (n[2], n[3]), (n[3], n[0])]
                        }
                    };

                    for (a, b) in ids {
                        let edge_key = if a < b { (a, b) } else { (b, a) };
                        // add face normal to edge
                        edge_normals.entry(edge_key).or_default().push(*normal);
                    }
                }
            }
        }

        // filter by angle
        let threshold_deg = 30.0_f64;
        let threshold_val = threshold_deg.to_radians().cos();

        self.wireframe_lines = edge_normals
            .iter()
            .filter_map(|((id_a, id_b), normals)| {
                match normals.len() {
                    1 => {
                        // Outer Edge, draw always (open mesh)
                        let n_a = self.nodes.get(id_a)?;
                        let n_b = self.nodes.get(id_b)?;
                        Some((Point3::from(n_a), Point3::from(n_b)))
                    }
                    2 => {
                        let n1 = normals[0];
                        let n2 = normals[1];
                        let dot = n1.dot(&n2);

                        if dot.abs() < threshold_val {
                            let n_a = self.nodes.get(id_a)?;
                            let n_b = self.nodes.get(id_b)?;
                            Some((Point3::from(n_a), Point3::from(n_b)))
                        } else {
                            None
                        }
                    }
                    _ => None, // Non-manifold, ignore
                }
            })
            .collect();
    }

    fn compute_face_normal(&self, face: &Face) -> Option<Unit<Vector3<f64>>> {
        let (i1, i2, i3) = match face {
            Face::Tri(idx) => (idx[0], idx[1], idx[2]),
            Face::Quad(idx) => (idx[0], idx[1], idx[2]), // Annahme: Planar, erste 3 Punkte reichen
        };

        let p1 = self.nodes.get(&i1)?;
        let p2 = self.nodes.get(&i2)?;
        let p3 = self.nodes.get(&i3)?;

        let v1 = Point3::new(p1.x, p1.y, p1.z);
        let v2 = Point3::new(p2.x, p2.y, p2.z);
        let v3 = Point3::new(p3.x, p3.y, p3.z);

        // Crossproduct: (p2 - p1) x (p3 - p1)
        let edge1 = v2 - v1;
        let edge2 = v3 - v1;

        Unit::try_new(edge1.cross(&edge2), 1e-6)
    }
}
