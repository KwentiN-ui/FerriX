use nalgebra::Point3;

use crate::solver::ids::NodeId;

#[derive(Debug, Clone)]
pub struct Node {
    pub id: NodeId,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Node {
    /// Attempts to parse a Node from the given stringslice
    pub fn parse_line(line: &str) -> Option<Self> {
        let delimited: Vec<&str> = line.split(',').map(str::trim).collect();
        if delimited.len() != 4 {
            return None;
        }
        let id = NodeId(delimited.first()?.parse::<usize>().ok()?);
        let x = delimited.get(1)?.parse::<f64>().ok()?;
        let y = delimited.get(2)?.parse::<f64>().ok()?;
        let z = delimited.get(3)?.parse::<f64>().ok()?;
        Some(Self { id, x, y, z })
    }
}

impl From<&Node> for Point3<f64> {
    fn from(n: &Node) -> Self {
        Point3::new(n.x, n.y, n.z)
    }
}
