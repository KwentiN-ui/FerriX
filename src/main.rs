use chrono::{Local, Utc};
use clap::Parser;
use rayon::ThreadPoolBuilder;

use crate::solver::{
    io::{vtk::VtkWriter, writer::ResultWriter},
    project::Project,
    results::StepResult,
    step::{static_step::StaticStep, steps::Step},
};

mod solver;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Build and install the global thread pool
    ThreadPoolBuilder::new()
        .num_threads(args.num_threads.unwrap_or(0)) // 0 means use default (num logical cores)
        .build_global()?;

    let project = Box::new(Project::from_jobname(&args.jobname, None)?);

    let mut all_results: Vec<StepResult> = Vec::new();

    println!("{}\n", Local::now().format("%d.%m.%y, %H:%M:%S"));
    println!("{LOGO}");
    println!("by Quentin Huss\n\n");

    println!("{}\n", project.get_info());

    let start_time = Utc::now();

    // Solver thread
    for (i, step_type) in project.steps.iter().enumerate() {
        let step_id = i + 1;
        match step_type {
            Step::StaticStep => {
                let mut step = StaticStep::new(project.clone());
                println!("--- Step {step_id}: StaticStep ---");
                match step.compute(step_id, &project.loads, &project.bcs) {
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
        let path = project
            .filepath
            .parent()
            .unwrap()
            .join(project.jobname() + ".vtk");

        // We use the mesh from the last step (assuming no remeshing)
        match writer.write(&path, &project.mesh.clone(), &all_results) {
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

    /// The number of threads to use for the analysis. Defaults to the number of logical cores.
    #[arg(short, long)]
    num_threads: Option<usize>,
}

const LOGO: &str = r" _____ _____ ____  ____  _ ___  _
/    //  __//  __\/  __\/ \\  \//
|  __\|  \  |  \/||  \/|| | \  /
| |   |  /_ |    /|    /| | /  \
\_/   \____\\_/\_\\_/\_\\_//__/\\
                                ";
