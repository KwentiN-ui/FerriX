use clap::Parser;

use crate::components::project::Project;

mod components;

fn main() {
    let args = Args::parse();

    let project = Project::from_jobname(&args.jobname, &args.preprocessed_output);
    match project {
        Ok(project) => {
            println!("Ok!");
        }
        Err(e) => eprintln!("{e}"),
    }
}

pub struct StaticStep;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
struct Args {
    /// The jobname. `CalculiX` will look for `<jobname>.inp`.
    jobname: String,

    /// Output path to write the preprocessed .inp file into. Useful for debugging
    #[arg(short, long)]
    preprocessed_output: Option<String>,
}
