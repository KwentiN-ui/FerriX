use clap::Parser;
use crossterm::event;
use std::{sync::mpsc, thread, time::Duration};

use crate::{
    solver::{
        project::Project,
        step::{static_step::StaticStep, steps::Step},
    },
    tui::{
        app::{App, AppEvent},
        setup, ui,
    },
};

mod solver;
mod tui;

#[allow(clippy::collapsible_if)]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let mut app = App::new(Project::from_jobname(&args.jobname, None)?);

    // Channel & Worker Thread
    let (tx, rx) = mpsc::channel();
    let tx_solver = tx.clone();

    let steps = app.project.steps.clone();
    // clone is fine here, it's just the pointer
    let input = app.project.input.clone();
    let mesh = app.project.mesh.clone();
    thread::spawn(move || {
        // Solver thread
        for (i, step_type) in steps.iter().enumerate() {
            match step_type {
                Step::StaticStep => {
                    let mut step = StaticStep::new(input.clone(), mesh.clone());
                    let _ = tx_solver
                        .send(AppEvent::SolverLog(format!("--- Step {i}: StaticStep ---")));
                    step.compute(&tx_solver);
                }
            }
        }
        tx_solver.send(AppEvent::AnalysisFinished).unwrap();
    });

    // TUI
    let mut terminal = setup::init()?;

    // Main Loop
    while !app.should_quit {
        terminal.draw(|f| ui::draw(f, &app))?;

        // Event handling
        if event::poll(Duration::from_millis(50))? {
            if let event::Event::Key(key) = event::read()? {
                app.update(AppEvent::Input(key));
            }
        }

        // Check for messages
        while let Ok(msg) = rx.try_recv() {
            app.update(msg);
        }
    }

    // cleanup
    setup::restore()?;
    Ok(())
}

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
struct Args {
    /// The jobname. `CalculiX` will look for `<jobname>.inp`.
    jobname: String,

    /// Output path to write the preprocessed .inp file into. Useful for debugging
    #[arg(short, long)]
    preprocessed_output: Option<String>,
}
