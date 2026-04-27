//! `FerriX` CLI
//!
//! The main entry point for the `FerriX` FEA solver. This executable handles command-line arguments,
//! initializes the simulation environment, and executes the solver steps.

use std::sync::Arc;

use chrono::{Local, Utc};
use clap::Parser;
use ferrix::solver::{io::OutputFormat, project::Project, state::SolutionState, time::SolverTime};
use rayon::ThreadPoolBuilder;

/// Default number of threads used for parallel computations.
const DEFAULT_THREADS: usize = 4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Build and install the global thread pool
    ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .build_global()?;

    let project = Arc::new(Project::from_jobname(&args.jobname, None)?);

    println!("{}\n", Local::now().format("%d.%m.%y, %H:%M:%S"));
    println!("{LOGO}");
    println!("by Quentin Huss\n\n");

    println!("{}", project.get_info());
    println!(
        "Using {} Thread(s)\nTry --help for more info\n\n",
        rayon::current_num_threads()
    );

    let start_time = Utc::now();
    let mut simulation_time = SolverTime::new();

    // --- Main solver loop ---
    let num_nodes = project.mesh.nodes.len();
    let num_dofs = num_nodes * 3;
    let mut solution_state = SolutionState::new(num_dofs, num_nodes);
    solution_state.initialize(&project);

    let writer = args.output_format.get_writer(project.clone());

    writer.init()?;

    for (i, step) in project.steps.iter().enumerate() {
        // Each step needs to call the simulation_time methods internally to advance the time properly!
        let step_id = i + 1;
        if let Err(e) = step.solve(
            step_id,
            &project,
            &mut solution_state,
            &*writer,
            &mut simulation_time,
        ) {
            eprintln!("Error occurred in step {step_id}: {e}\n\n");
            break;
        }
        println!("Step {step_id} completed.\n\n");
    }
    // --- End main solver loop ---

    let _ = writer.finish();

    println!("Job finished");
    println!(
        "________________________________________\nTotal FerriX Time: {:.3}s\n________________________________________",
        (Utc::now() - start_time).as_seconds_f64()
    );
    Ok(())
}

/// Command-line arguments for the `FerriX` solver.
#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
pub struct Args {
    /// The jobname. `FerriX` will look for `<jobname>.inp`.
    jobname: String,

    /// The number of threads to use for the analysis. 0 will use system CPU count.
    #[arg(short, long, default_value_t = DEFAULT_THREADS)]
    num_threads: usize,

    /// The format in which results will be saved.
    #[arg(short, long, default_value_t = OutputFormat::Vtk)]
    output_format: OutputFormat,

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
