//! Node definitions and parsing.
//!
//! A node represents a single point in 3D space with a unique identifier.

use crate::solver::ids::NodeId;
use nalgebra::Point3;

/// A single node in the finite element mesh.
#[derive(Debug, Clone)]
pub struct Node {
    /// Unique identifier for the node.
    pub id: NodeId,
    /// X-coordinate.
    pub x: f64,
    /// Y-coordinate.
    pub y: f64,
    /// Z-coordinate.
    pub z: f64,
}

impl Node {
    /// Parses a node from a comma-separated string line (e.g., "ID, X, Y, Z").
    ///
    /// Returns `None` if the line does not contain exactly 4 valid numeric fields.
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
