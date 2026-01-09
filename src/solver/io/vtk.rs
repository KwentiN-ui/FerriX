use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::writer::ResultWriter;
use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::results::{FieldType, StepResult};
use nalgebra::{Matrix3, SymmetricEigen};

pub struct VtkWriter;

impl ResultWriter for VtkWriter {
    #[allow(clippy::too_many_lines)]
    fn write(
        &self,
        path: &Path,
        mesh: &Mesh,
        results: &[StepResult],
        nodal_output: &[String],
        element_output: &[String],
    ) -> Result<(), Box<dyn Error>> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        // --- Header ---
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "CCX-RS Analysis Results")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET UNSTRUCTURED_GRID")?;

        // --- POINTS (Nodes) ---
        let num_nodes = mesh.index_to_node_id.len();
        writeln!(w, "POINTS {num_nodes} float")?;

        for &node_id in &mesh.index_to_node_id {
            let node = &mesh.nodes[&node_id];
            writeln!(w, "{} {} {}", node.x, node.y, node.z)?;
        }

        // --- CELLS (Elements) ---
        let mut cell_data: Vec<String> = Vec::new();
        let mut cell_types: Vec<u8> = Vec::new();
        let mut list_size = 0;

        for elem in mesh.elements.values() {
            match elem {
                Element::C3D4(_, nodes) => {
                    let idx0 = mesh.get_index_for_node_id(nodes[0]).unwrap();
                    let idx1 = mesh.get_index_for_node_id(nodes[1]).unwrap();
                    let idx2 = mesh.get_index_for_node_id(nodes[2]).unwrap();
                    let idx3 = mesh.get_index_for_node_id(nodes[3]).unwrap();
                    cell_data.push(format!("4 {idx0} {idx1} {idx2} {idx3}"));
                    cell_types.push(10); // VTK_TETRA
                    list_size += 5;
                }
                Element::C3D20(_, _nodes) => todo!("This is not supported yet!"),
            }
        }

        writeln!(w, "CELLS {} {}", cell_data.len(), list_size)?;
        for line in &cell_data {
            writeln!(w, "{line}")?;
        }

        writeln!(w, "CELL_TYPES {}", cell_types.len())?;
        for t in &cell_types {
            writeln!(w, "{t}")?;
        }

        // --- POINT_DATA ---
        if !results.is_empty() && (!nodal_output.is_empty() || !element_output.is_empty()) {
            writeln!(w, "POINT_DATA {num_nodes}")?;

            for step in results {
                // --- Displacement ---
                if nodal_output.contains(&"U".to_string()) {
                    if let Some(field) = step.nodal_results.iter().find(|f| f.field_type == FieldType::Displacement) {
                        writeln!(w, "VECTORS U float")?;
                        for &node_id in &mesh.index_to_node_id {
                            if let Some(val) = field.data.get(&node_id) {
                                writeln!(w, "{} {} {}", val[0], val[1], val[2])?;
                            } else {
                                writeln!(w, "0.0 0.0 0.0")?;
                            }
                        }
                    }
                }
                
                // --- Stress ---
                if element_output.contains(&"S".to_string()) {
                    if let Some(field) = step.nodal_results.iter().find(|f| f.field_type == FieldType::Stress) {
                        writeln!(w, "SCALARS S float 6")?;
                        writeln!(w, "LOOKUP_TABLE default")?;
                        for &node_id in &mesh.index_to_node_id {
                            if let Some(val) = field.data.get(&node_id) {
                                writeln!(w, "{} {} {} {} {} {}", val[0], val[1], val[2], val[3], val[4], val[5])?;
                            } else {
                                writeln!(w, "0.0 0.0 0.0 0.0 0.0 0.0")?;
                            }
                        }
                    }
                }

                // --- Strain ---
                if element_output.contains(&"E".to_string()) {
                    if let Some(field) = step.nodal_results.iter().find(|f| f.field_type == FieldType::Strain) {
                        writeln!(w, "SCALARS E float 6")?;
                        writeln!(w, "LOOKUP_TABLE default")?;
                        for &node_id in &mesh.index_to_node_id {
                            if let Some(val) = field.data.get(&node_id) {
                                writeln!(w, "{} {} {} {} {} {}", val[0], val[1], val[2], val[3], val[4], val[5])?;
                            } else {
                                writeln!(w, "0.0 0.0 0.0 0.0 0.0 0.0")?;
                            }
                        }
                    }
                }
                
                // --- MISES & TRESCA ---
                if element_output.contains(&"S".to_string()) {
                    let stress_field = step
                        .nodal_results
                        .iter()
                        .find(|f| f.field_type == FieldType::Stress);

                    if let Some(stress_field) = stress_field {
                        writeln!(w, "SCALARS MISES float 1")?;
                        writeln!(w, "LOOKUP_TABLE default")?;
                        for &node_id in &mesh.index_to_node_id {
                            if let Some(val) = stress_field.data.get(&node_id) {
                                let s11 = val[0];
                                let s22 = val[1];
                                let s33 = val[2];
                                let s12 = val[3];
                                let s23 = val[4];
                                let s13 = val[5];
                                let mises = ((s11 - s22).powi(2)
                                    + (s22 - s33).powi(2)
                                    + (s33 - s11).powi(2)
                                    + 6.0 * (s12.powi(2) + s23.powi(2) + s13.powi(2)))
                                    / 2.0;
                                writeln!(w, "{}", mises.sqrt())?;
                            } else {
                                writeln!(w, "0.0")?;
                            }
                        }

                        writeln!(w, "SCALARS TRESCA float 1")?;
                        writeln!(w, "LOOKUP_TABLE default")?;
                        for &node_id in &mesh.index_to_node_id {
                            if let Some(val) = stress_field.data.get(&node_id) {
                                let s11 = val[0];
                                let s22 = val[1];
                                let s33 = val[2];
                                let s12 = val[3];
                                let s23 = val[4];
                                let s13 = val[5];
                                let stress_matrix = Matrix3::new(
                                    s11, s12, s13, s12, s22, s23, s13, s23, s33,
                                );
                                let eigenvalues = SymmetricEigen::new(stress_matrix).eigenvalues;
                                let max_eigenvalue = eigenvalues.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                                let min_eigenvalue = eigenvalues.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                                let tresca = max_eigenvalue - min_eigenvalue;
                                writeln!(w, "{tresca}")?;
                            } else {
                                writeln!(w, "0.0")?;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
