use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::writer::ResultWriter;
use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::mesh_lib::mesh::Mesh;
use crate::solver::results::{FieldType, StepResult};

pub struct FrdWriter;

impl ResultWriter for FrdWriter {
    fn write(
        &self,
        path: &Path,
        mesh: &Mesh,
        results: &[StepResult],
    ) -> Result<(), Box<dyn Error>> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        // --- HEADER BLOCK (1C & 1U) ---
        // 1C marks the start. PrePoMax needs this.
        writeln!(w, "    1C")?;
        writeln!(w, "    1UTITLE        CCX-RS Analysis Result")?;
        writeln!(w, "    1UPGM          CCX-RS")?;

        // --- Block 1P: Nodes ---
        // Format: 1P, num_nodes, 1 (Format Code 1 implies: I2, I10, 3E12.5)
        writeln!(w, "    1P,{:10},    1", mesh.nodes.len())?;

        for node in mesh.nodes.values() {
            // Fix: Added " -1" (I2) at the start
            // Format: Key(-1), ID(I10), X, Y, Z (E12.5)
            writeln!(
                w,
                " -1{:10}{:12.5E}{:12.5E}{:12.5E}",
                node.id, node.x, node.y, node.z
            )?;
        }

        // --- Block 2C: Elements ---
        // Format: 2C, num_elems, 1 (Format Code 1 implies: I2, I10, I5, I5, I5, Nodes...)
        writeln!(w, "    2C,{:10},    1", mesh.elements.len())?;

        for elem in mesh.elements.values() {
            let (id, type_code, nodes) = match elem {
                // 9 = Tetra4 (C3D4)
                Element::C3D4(id, n) => (*id, 9, n.as_slice()),
                _ => continue,
            };

            // Fix: Added " -1" at the start
            // Format: Key(-1), ID(I10), Type(I5), Region(I5), 0(I5), Nodes(I10)...
            write!(w, " -1{:10}{:5}{:5}    0", id, type_code, 1)?;

            for &n_id in nodes {
                write!(w, "{:10}", n_id)?;
            }
            writeln!(w)?;
        }

        // --- Block 100: Nodal Results ---
        for step in results {
            for field in &step.nodal_results {
                if field.field_type == FieldType::Displacement {
                    writeln!(w, "    100CL")?;
                    writeln!(w, "{}_Step{}", field.name, step.step_id)?;
                    writeln!(
                        w,
                        " {:12.5E}, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00",
                        step.time_increment
                    )?;

                    writeln!(w, "TYPE")?;
                    writeln!(w, " DISPLACEMENT")?;

                    // Header Definition
                    writeln!(w, "    1    DISP    1    1    3    1")?;
                    writeln!(w, "    4  D1  D2  D3")?;

                    // Data Lines
                    // Format: Key(-1), ID(I10), Values(E12.5)
                    for (node_id, val) in &field.data {
                        writeln!(
                            w,
                            " -1{:10}{:12.5E}{:12.5E}{:12.5E}",
                            node_id, val[0], val[1], val[2]
                        )?;
                    }
                    writeln!(w, "    -1")?; // End Dataset Marker
                }
            }
        }

        writeln!(w, "    9999")?; // Global End of File
        Ok(())
    }
}
