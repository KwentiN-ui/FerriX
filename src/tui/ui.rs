use nalgebra::Point3;
use ratatui::{
    prelude::*,
    symbols::Marker,
    widgets::{
        Block, Borders, List, ListItem,
        canvas::{Canvas, Context, Line, Points},
    },
};

use crate::tui::app::App;

pub fn draw(f: &mut Frame, app: &App) {
    let chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(50), Constraint::Min(50)])
        .split(f.area());

    // Mesh canvas
    let jobname = app.project.jobname();
    let canvas = Canvas::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(format!(" Mesh Preview - {jobname} "))
                .title_bottom(" Zoom - X/Y, Cycle Up - U, Rotate - Up/Down/Left/Right "),
        )
        .marker(Marker::Braille) // Wichtig für hohe Auflösung!
        .x_bounds([-1.0, 1.0])
        .y_bounds([-1.0, 1.0])
        .paint(|ctx| {
            draw_mesh(ctx, app);
        });

    f.render_widget(canvas, chunks[0]);

    // Log
    let log_items: Vec<ListItem> = app
        .logs
        .iter()
        .map(|msg| ListItem::new(msg.as_str()))
        .collect();

    let logs = List::new(log_items)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Solver Logs "),
        )
        .highlight_style(Style::default().add_modifier(Modifier::BOLD));

    f.render_widget(logs, chunks[1]);
}

fn draw_mesh(ctx: &mut Context, app: &App) {
    let project = &app.project;
    let mvp = app.camera.build_view_projection_matrix();

    for (p_start, p_end) in &project.mesh.wireframe_lines {
        let s_ndc = mvp.transform_point(p_start);
        let e_ndc = mvp.transform_point(p_end);

        // Simple Culling Check
        if s_ndc.z.abs() <= 1.0 && e_ndc.z.abs() <= 1.0 {
            ctx.draw(&Line {
                x1: s_ndc.x,
                y1: s_ndc.y,
                x2: e_ndc.x,
                y2: e_ndc.y,
                color: Color::Yellow,
            });
        }
    }
}
