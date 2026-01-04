use crate::components::mesh_lib::{elements::element::Element, mesh::Mesh};

pub struct C3D20;

impl Element for C3D20 {
    fn name() -> String {
        "C3D20".to_string()
    }
}
