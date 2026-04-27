//! VTK result exporter.
//!
//! This module implements the `ResultWriter` trait for the VTK format,
//! allowing results to be visualized in software like `ParaView`.

use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use super::writer::ResultWriter;
use crate::solver::mesh_lib::elements::element::Element;
use crate::solver::project::Project;
use crate::solver::results::{FieldType, IncResult};
use crate::solver::time::SolverTime;
use nalgebra::{Matrix3, SymmetricEigen};

/// A result writer that exports simulation data to VTK files.
///
/// It generates unstructured grid (.vtk) files for each increment and a
/// ParaView-compatible (.vtk.series) file to manage the time steps.
pub struct VtkWriter {
    project: Arc<Project>,
}

impl VtkWriter {
    /// Creates a new `VtkWriter`.
    #[must_use]
    pub fn new(project: Arc<Project>) -> Self {
        Self { project }
    }
    fn dirpath(&self) -> PathBuf {
        let jobname = self
            .project
            .jobname()
            .unwrap_or_else(|_| "Unknown".to_string());
        self.project
            .job_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(jobname)
    }
}

const SERIES_FILENAME: &str = "results.vtk.series";

impl ResultWriter for VtkWriter {
    /// Is called at the beginning of analysis. Can be used to setup directories etc.
    ///
    /// # Errors
    /// Returns an error if the initialization fails (e.g. directory creation).
    fn init(&self) -> Result<(), Box<dyn Error>> {
        let dirpath = self.dirpath();
        fs::create_dir_all(&dirpath)?;

        let series = File::create(dirpath.join(SERIES_FILENAME))?;
        let mut w = BufWriter::new(series);
        writeln!(
            w,
            "{{
          \"file-series-version\" : \"1.0\",
          \"files\" : ["
        )?;
        Ok(())
    }

    /// Writes the results of an increment.
    ///
    /// # Errors
    /// Returns an error if the writing fails.
    #[allow(clippy::too_many_lines)]
    fn write_increment(
        &self,
        inc_result: &IncResult,
        timer: &SolverTime,
    ) -> Result<(), Box<dyn Error>> {
        // Create output directory

        // Creates a new folder with the jobname and writes increments into it
        let filename = format!["step_{}_{}.vtk", inc_result.step_id, inc_result.inc_id];
        let file = File::create(self.dirpath().join(&filename))?;
        let mut w = BufWriter::new(file);

        // --- Header ---
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "CCX-RS Analysis Results")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET UNSTRUCTURED_GRID")?;

        // --- POINTS (Nodes) ---
        let mesh = self.project.mesh.clone();
        let num_nodes = mesh.index_to_node_id.len();
        writeln!(w, "POINTS {num_nodes} float")?;

        for &node_id in &mesh.index_to_node_id {
            if let Some(node) = mesh.nodes.get(&node_id) {
                writeln!(w, "{} {} {}", node.x, node.y, node.z)?;
            } else {
                writeln!(w, "0.0 0.0 0.0")?;
            }
        }

        // --- CELLS (Elements) ---
        let mut cell_data: Vec<String> = Vec::new();
        let mut cell_types: Vec<u8> = Vec::new();
        let mut list_size = 0;

        for elem in mesh.elements.values() {
            match elem {
                Element::C3D4(_, nodes) => {
                    if let (Some(idx0), Some(idx1), Some(idx2), Some(idx3)) = (
                        mesh.get_index_for_node_id(nodes[0]),
                        mesh.get_index_for_node_id(nodes[1]),
                        mesh.get_index_for_node_id(nodes[2]),
                        mesh.get_index_for_node_id(nodes[3]),
                    ) {
                        cell_data.push(format!("4 {idx0} {idx1} {idx2} {idx3}"));
                        cell_types.push(10); // VTK_TETRA
                        list_size += 5;
                    }
                }
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
        writeln!(w, "POINT_DATA {num_nodes}")?;

        // --- Displacement ---
        if let Some(field) = inc_result
            .nodal_results
            .iter()
            .find(|f| f.field_type == FieldType::Displacement)
        {
            writeln!(w, "VECTORS U float")?;
            for &node_id in &mesh.index_to_node_id {
                if let Some(val) = field.data.get(&node_id) {
                    writeln!(w, "{} {} {}", val[0], val[1], val[2])?;
                } else {
                    writeln!(w, "0.0 0.0 0.0")?;
                }
            }
        }

        // --- Stress ---
        if let Some(field) = inc_result
            .nodal_results
            .iter()
            .find(|f| f.field_type == FieldType::Stress)
        {
            writeln!(w, "SCALARS S float 6")?;
            writeln!(w, "LOOKUP_TABLE default")?;
            for &node_id in &mesh.index_to_node_id {
                if let Some(val) = field.data.get(&node_id) {
                    writeln!(
                        w,
                        "{} {} {} {} {} {}",
                        val[0], val[1], val[2], val[3], val[4], val[5]
                    )?;
                } else {
                    writeln!(w, "0.0 0.0 0.0 0.0 0.0 0.0")?;
                }
            }
        }

        // --- Strain ---
        if let Some(field) = inc_result
            .nodal_results
            .iter()
            .find(|f| f.field_type == FieldType::Strain)
        {
            writeln!(w, "SCALARS E float 6")?;
            writeln!(w, "LOOKUP_TABLE default")?;
            for &node_id in &mesh.index_to_node_id {
                if let Some(val) = field.data.get(&node_id) {
                    writeln!(
                        w,
                        "{} {} {} {} {} {}",
                        val[0], val[1], val[2], val[3], val[4], val[5]
                    )?;
                } else {
                    writeln!(w, "0.0 0.0 0.0 0.0 0.0 0.0")?;
                }
            }
        }

        // --- Temperature ---
        if let Some(field) = inc_result
            .nodal_results
            .iter()
            .find(|f| f.field_type == FieldType::Temperature)
        {
            writeln!(w, "SCALARS NT float 1")?;
            writeln!(w, "LOOKUP_TABLE default")?;
            for &node_id in &mesh.index_to_node_id {
                if let Some(val) = field.data.get(&node_id) {
                    writeln!(w, "{}", val[0])?;
                } else {
                    writeln!(w, "0.0")?;
                }
            }
        }

        // --- MISES & TRESCA ---
        let stress_field = inc_result
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
                    let stress_matrix = Matrix3::new(s11, s12, s13, s12, s22, s23, s13, s23, s33);
                    let eigenvalues = SymmetricEigen::new(stress_matrix).eigenvalues;
                    let max_eigenvalue =
                        eigenvalues.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                    let min_eigenvalue = eigenvalues.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                    let tresca = max_eigenvalue - min_eigenvalue;
                    writeln!(w, "{tresca}")?;
                } else {
                    writeln!(w, "0.0")?;
                }
            }
        }

        // add file to series
        let series = OpenOptions::new()
            .append(true)
            .open(self.dirpath().join(SERIES_FILENAME))?;
        let mut w = BufWriter::new(series);
        writeln!(
            w,
            "{{ \"name\" : \"{}\", \"time\" : {} }},",
            &filename,
            timer.global_time()
        )?;

        Ok(())
    }

    /// Is called at the very end of the analysis. Can be used for cleanup, etc.
    ///
    /// # Errors
    /// Returns an error if the finish operation fails.
    fn finish(&self) -> Result<(), Box<dyn Error>> {
        let series = OpenOptions::new()
            .append(true)
            .open(self.dirpath().join(SERIES_FILENAME))?;

        let mut w = BufWriter::new(series);
        writeln!(w, "  ]\n}}")?; // Schließt das "files"-Array und das Hauptobjekt
        Ok(())
    }
}
