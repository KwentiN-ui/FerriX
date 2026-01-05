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

        // --- Block 1P: Nodes (Global, Static Topology) ---
        writeln!(w, "    1P,{},    1", mesh.nodes.len())?;
        for node in mesh.nodes.values() {
            writeln!(
                w,
                "{:10} {:12.5E} {:12.5E} {:12.5E}",
                node.id, node.x, node.y, node.z
            )?;
        }

        // --- Block 2C: Elements (Global, Static Topology) ---
        writeln!(w, "    2C,{},    1", mesh.elements.len())?;
        for elem in mesh.elements.values() {
            let (id, type_code, nodes) = match elem {
                Element::C3D4(id, n) => (*id, 9, n.as_slice()),
                _ => continue,
            };

            write!(w, "{:10}{:5}{:5}    0", id, type_code, 0)?;
            for &n_id in nodes {
                write!(w, "{n_id:10}")?;
            }
            writeln!(w)?;
        }

        // --- Block 100: Nodal Results (Per Step) ---
        for step in results {
            for field in &step.nodal_results {
                if field.field_type == FieldType::Displacement {
                    // 100CL indicates a generic dataset
                    writeln!(w, "    100CL")?;

                    // Dataset Name (displayed in cgx)
                    writeln!(w, "{}_Step{}", field.name, step.step_id)?;

                    // Value Block: Important for Time/Frequency mapping
                    // We place the time_increment here so cgx can order them
                    writeln!(
                        w,
                        " {:12.5E}, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00",
                        step.time_increment
                    )?;

                    writeln!(w, "TYPE")?;
                    writeln!(w, " DISPLACEMENT")?;

                    // Header: Key, Name, Menu, ICtype(1=Node), Comp(3), IRtype(1=Float)
                    writeln!(w, "    1    DISP    1    1    3    1")?;
                    writeln!(w, "    4  D1  D2  D3")?; // Component names

                    // Data Lines
                    for (node_id, val) in &field.data {
                        writeln!(
                            w,
                            "{:10}{:12.5E}{:12.5E}{:12.5E}",
                            node_id, val[0], val[1], val[2]
                        )?;
                    }
                    writeln!(w, "    -1")?; // End Dataset
                }
            }
        }

        writeln!(w, "    9999")?; // End of File
        Ok(())
    }
}
