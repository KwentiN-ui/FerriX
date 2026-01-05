use ndarray::{Array2, array};
use std::error::Error;
use std::str::FromStr;

use strum_macros::{EnumDiscriminants, EnumString};

/// <https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node194.html>
/// strum automatically generates a String-enum `ElementType` based on these definitions.
#[derive(EnumDiscriminants, Debug, Clone)]
#[strum_discriminants(derive(Hash, EnumString))]
#[strum_discriminants(name(ElementType))]
pub enum Element {
    // General 3D-Solids
    /// 4-node linear tetrahedral element
    C3D4(usize, [usize; 4]),
    /// 6-node linear triangular prism element
    C3D6(usize, [usize; 6]),
    /// 3D 20-node quadratic isoparametric element
    C3D20(usize, [usize; 20]),

    // Shell elements
    /// S8 (8-node quadratic shell element)
    S8(usize, [usize; 8]),
}

impl Element {
    pub fn parse_type_str_from_line(line: &str) -> Result<String, Box<dyn Error>> {
        Ok(line
            .split(',')
            .map(str::trim)
            .nth(1)
            .ok_or("Invalid element definition on line {line_nr}")?
            .split('=')
            .next_back()
            .ok_or("Invalid element definition on line {line_nr}")?
            .to_string())
    }
    /// Create an Element from a line. This function panics if it's not able to do so.
    pub fn parse_line(type_name: &str, line: &str) -> Self {
        let nums: Vec<usize> = line
            .split(',')
            .map(|s| s.trim().parse().expect("Integer conversion failed"))
            .collect();

        let (&id, nodes) = nums.split_first().expect("Line empty");

        // String -> ElementType
        let elem_type = ElementType::from_str(type_name)
            .unwrap_or_else(|_| panic!("Unknown element definition: {type_name}"));

        // Local Macro for array casting
        macro_rules! to_arr {
            ($n:expr) => {
                nodes
                    .try_into()
                    .expect(concat!("Wrong node count for ", stringify!($n)))
            };
        }

        match elem_type {
            ElementType::C3D4 => Element::C3D4(id, to_arr!(C3D4)),
            ElementType::C3D6 => Element::C3D6(id, to_arr!(C3D6)),
            ElementType::C3D20 => Element::C3D20(id, to_arr!(C3D20)),
            ElementType::S8 => Element::S8(id, to_arr!(S8)),
        }
    }

    pub fn get_faces(&self) -> Vec<Face> {
        match self {
            // C3D4: 4 Triangle Faces
            Element::C3D4(_, n) => vec![
                Face::tri(n[0], n[1], n[2]),
                Face::tri(n[0], n[1], n[3]),
                Face::tri(n[1], n[2], n[3]),
                Face::tri(n[2], n[0], n[3]),
            ],
            // C3D6: 2 Triangle Faces, Top and Bottom
            Element::C3D6(_, n) => vec![
                Face::tri(n[0], n[1], n[2]), // Unten
                Face::tri(n[3], n[4], n[5]), // Oben
                Face::quad(n[0], n[1], n[4], n[3]),
                Face::quad(n[1], n[2], n[5], n[4]),
                Face::quad(n[2], n[0], n[3], n[5]),
            ],
            // C3D20: 6 Quad Faces
            Element::C3D20(_, n) => vec![
                Face::quad(n[0], n[1], n[2], n[3]), // Unten
                Face::quad(n[4], n[5], n[6], n[7]), // Oben
                Face::quad(n[0], n[1], n[5], n[4]), // Seiten...
                Face::quad(n[1], n[2], n[6], n[5]),
                Face::quad(n[2], n[3], n[7], n[6]),
                Face::quad(n[3], n[0], n[4], n[7]),
            ],
            // S8: Is a single Face
            Element::S8(_, n) => vec![Face::quad(n[0], n[1], n[2], n[3])],
        }
    }

    pub fn get_id(&self) -> usize {
        match self {
            Element::C3D20(id, _)
            | Element::C3D6(id, _)
            | Element::C3D4(id, _)
            | Element::S8(id, _) => *id,
        }
    }

    /// Get global Node IDs
    pub fn get_node_ids(&self) -> &[usize] {
        match self {
            Element::C3D4(_, n) => n,
            Element::C3D6(_, n) => n,
            Element::C3D20(_, n) => n,
            Element::S8(_, n) => n,
        }
    }

    pub fn integration_points(&self) -> Vec<GaussPoint> {
        match self {
            Element::C3D4(..) => c3d4_gauss(),
            Element::C3D20(..) => c3d20_gauss(),
            _ => todo!("This is not implement yet..."),
        }
    }

    pub fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> (Vec<f64>, Array2<f64>) {
        match self {
            Element::C3D4(..) => shape_func_c3d4(xi, eta, zeta),
            Element::C3D20(..) => shape_func_c3d20(xi, eta, zeta),
            _ => todo!("Shape functions missing"),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GaussPoint {
    pub coords: [f64; 3], // xi, eta, zeta
    pub weight: f64,
}

use std::cmp::{max, min};

use crate::solver::mesh_lib::elements::{
    c3d4::{c3d4_gauss, shape_func_c3d4},
    c3d20::{c3d20_gauss, shape_func_c3d20},
};

// Helper for Drawing the Mesh
#[derive(Debug, Hash, Eq, PartialEq, Clone)]
pub enum Face {
    Tri([usize; 3]),
    Quad([usize; 4]),
}

impl Face {
    fn tri(a: usize, b: usize, c: usize) -> Self {
        let mut arr = [a, b, c];
        arr.sort_unstable(); // make face canonical
        Face::Tri(arr)
    }

    fn quad(a: usize, b: usize, c: usize, d: usize) -> Self {
        let mut arr = [a, b, c, d];
        arr.sort_unstable();
        Face::Quad(arr)
    }
}
