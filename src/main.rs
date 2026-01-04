use clap::Parser;

use crate::components::project::Project;

mod components;

fn main() {
    let args = Args::parse();

    let project = Project::from_jobname(&args.jobname);
    match project {
        Ok(project) => {
            println!("Ok!");
        }
        Err(e) => println!("{e}"),
    }
}

pub struct StaticStep;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
struct Args {
    /// The jobname. `CalculiX` will look for `<jobname>.inp`.
    jobname: String,
}
