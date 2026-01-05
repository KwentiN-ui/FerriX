use std::error::Error;

/// <https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node194.html>
#[derive(Debug)]
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
            .ok_or("Invalid Element definition on line {line_nr}")?
            .split('=')
            .next_back()
            .ok_or("Invalid Element definition on line {line_nr}")?
            .to_string())
    }
    /// Create an Element from a line. This function panics if it's not able to do so.
    pub fn parse_line(type_name: &str, line: &str) -> Self {
        let nums: Vec<usize> = line
            .split(',')
            .map(|s| {
                s.trim()
                    .parse()
                    .expect("Could not convert {s} to an integer!")
            })
            .collect();

        // split ID and Nodes
        let (&id, nodes) = nums
            .split_first()
            .expect("Empty lines not possible after preprocessing.");

        match type_name {
            "C3D4" => {
                let arr: [usize; 4] = nodes
                    .try_into()
                    .expect("Wrong amount of elements for C3D4 definition!");
                Element::C3D4(id, arr)
            }
            "C3D6" => {
                let arr: [usize; 6] = nodes
                    .try_into()
                    .expect("Wrong amount of elements for C3D6 definition!");
                Element::C3D6(id, arr)
            }
            "C3D20" => {
                let arr: [usize; 20] = nodes
                    .try_into()
                    .expect("Wrong amount of elements for C3D20 definition!");
                Element::C3D20(id, arr)
            }
            "S8" => {
                let arr: [usize; 8] = nodes
                    .try_into()
                    .expect("Wrong amount of elements for S8 definition!");
                Element::S8(id, arr)
            }
            _ => panic!("Unknown element definition: {type_name}"),
        }
    }
}
