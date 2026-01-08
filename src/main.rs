use clap::Parser;

use crate::solver::{
    io::{vtk::VtkWriter, writer::ResultWriter},
    project::Project,
    results::StepResult,
    step::{static_step::StaticStep, steps::Step},
};

mod solver;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let project = Project::from_jobname(&args.jobname, None)?;

    let steps = project.steps.clone();
    let filepath = project.filepath.clone();
    // clone is fine here, it's just the pointer
    let input = project.input.clone();
    let mesh = project.mesh.clone();

    let mut all_results: Vec<StepResult> = Vec::new();

    println!("{}", project.get_info());

    // Solver thread
    for (i, step_type) in steps.iter().enumerate() {
        match step_type {
            Step::StaticStep(line) => {
                let mut step = StaticStep::new(input.clone(), mesh.clone(), *line);
                println!("--- Step {i}: StaticStep ---");
                match step.compute() {
                    Ok(res) => {
                        all_results.push(res);
                        println!("Step {i} completed.");
                    }
                    Err(e) => {
                        eprintln!(
                            "Error occured in step {i}: {e}\n\nAttempting to write results..."
                        );
                        break;
                    }
                }
            }
        }
    }
    // write results
    if !all_results.is_empty() {
        let writer = VtkWriter;
        let path = filepath.parent().unwrap().join("results.vtk");

        // We use the mesh from the last step (assuming no remeshing)
        match writer.write(&path, &mesh.clone(), &all_results) {
            Ok(()) => {
                println!("Written results to disk!");
            }
            Err(e) => {
                eprintln!("{e}");
            }
        }
    }
    println!("Analysis done, have a nice day!");

    Ok(())
}

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
pub struct Args {
    /// The jobname. `CalculiX` will look for `<jobname>.inp`.
    jobname: String,

    /// Output path to write the preprocessed .inp file into. Useful for debugging
    #[arg(short, long)]
    preprocessed_output: Option<String>,
}
