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
    pub camera: Camera,
}

impl App {
    pub fn new(jobname: String, project: Project) -> Self {
        let center = project.mesh.get_center();
        Self {
            should_quit: false,
            jobname,
            logs: vec!["Ready.".into()],
            project: Some(project),
            camera: Camera::new(center),
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
            KeyCode::Char('y') => self.camera.move_back(),
            KeyCode::Char('x') => self.camera.move_in(),
            KeyCode::Char('u') => self.camera.cycle_up(),
            KeyCode::Up => self.camera.rotate_up(),
            KeyCode::Down => self.camera.rotate_down(),
            KeyCode::Left => self.camera.rotate_left(),
            KeyCode::Right => self.camera.rotate_right(),
            _ => {}
        }
    }
}

use nalgebra::{Isometry3, Matrix4, Perspective3, Point3, Rotation3, Unit, Vector3};
#[derive(Debug, Default)]
pub struct Camera {
    pub pos: Point3<f64>,
    pub target: Point3<f64>,
    pub up: Vector3<f64>,
    world_up: Vector3<f64>,
    pub aspect: f64,
    /// FOV in rad
    pub fov: f64,
}

impl Camera {
    pub fn new(center: Point3<f64>) -> Self {
        Self {
            pos: Point3::new(20.0, 20.0, 20.0),
            target: center,
            up: Vector3::z(),
            world_up: Vector3::z(),
            aspect: 1.6,
            fov: std::f64::consts::PI / 4.0,
        }
    }

    pub fn move_back(&mut self) {
        self.pos -= (self.target - self.pos) * 0.5;
    }
    pub fn move_in(&mut self) {
        self.pos += (self.target - self.pos) * 0.5;
    }

    pub fn rotate_up(&mut self) {
        self.orbit_vertical(0.1);
    }

    pub fn rotate_down(&mut self) {
        self.orbit_vertical(-0.1);
    }

    pub fn rotate_left(&mut self) {
        self.orbit_horizontal(-0.1);
    }

    pub fn rotate_right(&mut self) {
        self.orbit_horizontal(0.1);
    }

    /// Rotates position AND up-vector around the local right-axis.
    fn orbit_vertical(&mut self, theta: f64) {
        let radius_vec = self.pos - self.target;
        let axis = radius_vec.cross(&self.up);

        if let Some(unit_axis) = Unit::try_new(axis, 1e-6) {
            let rotation = Rotation3::from_axis_angle(&unit_axis, -theta);
            self.pos = self.target + rotation * radius_vec;
            self.up = rotation * self.up;
        }
    }

    /// Rotates position around the GLOBAL world-up vector.
    fn orbit_horizontal(&mut self, theta: f64) {
        let radius_vec = self.pos - self.target;

        // Use world_up as axis, NOT self.up
        if let Some(unit_axis) = Unit::try_new(self.world_up, 1e-6) {
            let rotation = Rotation3::from_axis_angle(&unit_axis, theta);

            // Rotate position relative to target
            self.pos = self.target + rotation * radius_vec;

            // We MUST also rotate the camera's local up vector.
            // Otherwise, the camera twists relative to the turn direction.
            self.up = rotation * self.up;
        }
    }

    pub fn cycle_up(&mut self) {
        self.world_up = Vector3::new(self.world_up.z, self.world_up.x, self.world_up.y);
        self.up = self.world_up;
    }

    pub fn build_view_projection_matrix(&self) -> Matrix4<f64> {
        // 1. View Matrix (World -> Camera)
        let view = Isometry3::look_at_rh(&self.pos, &self.target, &self.up);

        let projection = Perspective3::new(self.aspect, self.fov, 0.1, 1000.0);

        projection.as_matrix() * view.to_homogeneous()
    }
}
