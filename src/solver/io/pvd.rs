use std::error::Error;
use std::path::PathBuf;
use crate::solver::io::vtk::VtkWriter;
use crate::solver::io::writer::ResultWriter;
use crate::solver::project::Project;
use crate::solver::results::StepResult;

#[derive(Default)]
pub struct PvdWriter {
    jobname: String,
    output_dir: PathBuf,
    vtk_writer: VtkWriter,
    pvd_content: String,
    project: Option<Project>,
    file_counter: usize,
}

impl PvdWriter {
    pub fn new() -> Self {
        Self::default()
    }
}

impl ResultWriter for PvdWriter {
    fn init(&mut self, project: &Project) -> Result<(), Box<dyn Error>> {
        self.project = Some(project.clone());
        self.jobname = project.jobname().to_string();
        self.output_dir = project.filepath.parent().unwrap().to_path_buf();
        self.pvd_content
            .push_str("<VTKFile type=\"Collection\" version=\"0.1\">\\n");
        self.pvd_content.push_str("  <Collection>\\n");
        Ok(())
    }

    fn write(&mut self, result: &StepResult) -> Result<(), Box<dyn Error>> {
        let project = self.project.as_ref().unwrap();
        let inc_filename = format!("{}_inc_{}.vtk", self.jobname, self.file_counter);
        let inc_path = self.output_dir.join(&inc_filename);

        // Use the VtkWriter to write the increment file
        self.vtk_writer.write_single_increment(
            &inc_path,
            &project.mesh,
            result,
            &project.nodal_output,
            &project.element_output,
        )?;

        // Add entry to PVD file content
        self.pvd_content.push_str(&format!(
            "    <DataSet timestep=\"{}\" group=\"\" part=\"0\" file=\"{}\"/>\\n",
            result.time_increment, inc_filename
        ));

        self.file_counter += 1;

        Ok(())
    }

    fn finish(&mut self) -> Result<(), Box<dyn Error>> {
        self.pvd_content.push_str("  </Collection>\\n");
        self.pvd_content.push_str("</VTKFile>\\n");

        let pvd_path = self.output_dir.join(format!("{}.pvd", self.jobname));
        std::fs::write(&pvd_path, &self.pvd_content)?;
        println!("Written PVD file to {}", pvd_path.display());
        Ok(())
    }
}

