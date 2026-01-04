use clap::Parser;

use crate::components::project::Project;

mod components;

fn main() {
    let args = Args::parse();

    let _project = Project::from_filepath(&args.inp_file);
}

pub struct StaticStep;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
struct Args {
    /// Filepath to the .inp file. The suffix is optional.
    inp_file: String,
}
