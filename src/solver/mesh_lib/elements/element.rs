use std::error::Error;
use std::str::FromStr;

use strum_macros::{EnumDiscriminants, EnumString};

/// <https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node194.html>
/// strum automatically generates a String-enum based on these definitions.
#[derive(EnumDiscriminants, Debug)]
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
}
