use chrono::{Local, Utc};
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

    println!("{}\n", Local::now().format("%d.%m.%y, %H:%M:%S"));
    println!("{LOGO}");
    println!("by Quentin Huss\n\n");

    println!("{}\n", project.get_info());

    let start_time = Utc::now();

    // Solver thread
    for (i, step_type) in steps.iter().enumerate() {
        let step_id = i + 1;
        match step_type {
            Step::StaticStep(line) => {
                let mut step = StaticStep::new(input.clone(), mesh.clone(), *line);
                println!("--- Step {step_id}: StaticStep ---");
                match step.compute(step_id) {
                    Ok(res) => {
                        all_results.push(res);
                        println!("Step {step_id} completed.\n\n");
                    }
                    Err(e) => {
                        eprintln!(
                            "Error occured in step {step_id}: {e}\n\nAttempting to write results..."
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
        let path = filepath.parent().unwrap().join(project.jobname() + ".vtk");

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
    println!("{:.2}s", (Utc::now() - start_time).as_seconds_f64());
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

const LOGO: &str = r" _____ _____ ____  ____  _ ___  _
/    //  __//  __\/  __\/ \\  \//
|  __\|  \  |  \/||  \/|| | \  /
| |   |  /_ |    /|    /| | /  \
\_/   \____\\_/\_\\_/\_\\_//__/\\
                                ";
