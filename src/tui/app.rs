use crossterm::event::{self, KeyCode, KeyEvent, KeyEventKind};

use crate::solver::project::Project;

pub enum AppEvent {
    Input(KeyEvent),
    SolverLog(String),
    SolverFinished,
}

pub struct App {
    pub should_quit: bool,
    pub jobname: String,
    pub logs: Vec<String>,
    pub project: Option<Project>,
}

impl App {
    pub fn new(jobname: String, project: Project) -> Self {
        Self {
            should_quit: false,
            jobname,
            logs: vec!["Ready.".into()],
            project: Some(project),
        }
    }

    pub fn update(&mut self, event: AppEvent) {
        match event {
            AppEvent::Input(key) => self.handle_key(key),
            AppEvent::SolverLog(msg) => self.logs.push(msg),
            AppEvent::SolverFinished => self.logs.push("Done.".into()),
        }
    }

    fn handle_key(&mut self, key: KeyEvent) {
        if key.kind != KeyEventKind::Press {
            return;
        }
        match key.code {
            KeyCode::Char('q') | KeyCode::Esc => self.should_quit = true,
            _ => {}
        }
    }
}
