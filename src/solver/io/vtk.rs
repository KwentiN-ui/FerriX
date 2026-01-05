use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::writer::ResultWriter;
use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::results::{FieldType, StepResult};

pub struct VtkWriter;

impl ResultWriter for VtkWriter {
    fn write(
        &self,
        path: &Path,
        mesh: &Mesh,
        results: &[StepResult],
    ) -> Result<(), Box<dyn Error>> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        // --- Header ---
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "CCX-RS Analysis Results")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET UNSTRUCTURED_GRID")?;

        // --- POINTS (Knoten) ---
        // VTK referenziert Knoten über ihren Index (0, 1, 2...) in dieser Liste.
        // Wir MÜSSEN daher zwingend unser `index_to_node_id` Mapping nutzen,
        // um die Reihenfolge konsistent zu halten.
        let num_nodes = mesh.index_to_node_id.len();
        writeln!(w, "POINTS {} float", num_nodes)?;

        for &node_id in &mesh.index_to_node_id {
            // Unsafe unwrap ist hier ok, da das Mapping konsistent sein muss
            let node = &mesh.nodes[&node_id];
            writeln!(w, "{} {} {}", node.x, node.y, node.z)?;
        }

        // --- CELLS (Elemente) ---
        // Wir müssen erst die Daten sammeln, um die Gesamtgröße für den Header zu kennen.
        // Format pro Zeile: "Anzahl_Knoten Index1 Index2 ..."
        // Indizes müssen die 0-basierten VTK-Indizes sein!

        let mut cell_data: Vec<String> = Vec::new();
        let mut cell_types: Vec<u8> = Vec::new();
        let mut list_size = 0; // Anzahl Einträge (Anzahl Elemente + Summe aller Knotenreferenzen)

        for elem in mesh.elements.values() {
            match elem {
                Element::C3D4(_, nodes) => {
                    // Mapping: NodeID -> VTK Index
                    // Wir nutzen das Mesh-Lookup
                    let idx0 = mesh.get_index_for_node_id(nodes[0]).unwrap();
                    let idx1 = mesh.get_index_for_node_id(nodes[1]).unwrap();
                    let idx2 = mesh.get_index_for_node_id(nodes[2]).unwrap();
                    let idx3 = mesh.get_index_for_node_id(nodes[3]).unwrap();

                    // 4 Knoten + 1 Count-Integer = 5 Einträge
                    cell_data.push(format!("4 {} {} {} {}", idx0, idx1, idx2, idx3));
                    cell_types.push(10); // VTK_TETRA
                    list_size += 5;
                }
                // Hier später C3D10, Hex8 etc. ergänzen
                _ => {}
            }
        }

        writeln!(w, "CELLS {} {}", cell_data.len(), list_size)?;
        for line in &cell_data {
            writeln!(w, "{}", line)?;
        }

        writeln!(w, "CELL_TYPES {}", cell_types.len())?;
        for t in &cell_types {
            writeln!(w, "{}", t)?;
        }

        // --- POINT_DATA (Ergebnisse) ---
        // Wenn wir Ergebnisse haben, schreiben wir sie
        if !results.is_empty() {
            writeln!(w, "POINT_DATA {}", num_nodes)?;

            for step in results {
                for field in &step.nodal_results {
                    if field.field_type == FieldType::Displacement {
                        // Wir benennen das Feld inklusive Step-Nummer
                        let field_name = format!("{}_Step{}", field.name, step.step_id);

                        writeln!(w, "VECTORS {} float", field_name)?;

                        // WICHTIG: Die Reihenfolge muss wieder exakt den POINTS entsprechen!
                        // Wir iterieren über `index_to_node_id` (0..N) und suchen den Wert.
                        for &node_id in &mesh.index_to_node_id {
                            if let Some(val) = field.data.get(&node_id) {
                                writeln!(w, "{} {} {}", val[0], val[1], val[2])?;
                            } else {
                                // Fallback, falls ein Knoten keine Ergebnisse hat (sollte nicht passieren)
                                writeln!(w, "0.0 0.0 0.0")?;
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }
}
