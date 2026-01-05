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
        .constraints([Constraint::Percentage(70), Constraint::Min(30)])
        .split(f.area());

    // Mesh canvas
    let canvas = Canvas::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Mesh Preview "),
        )
        .marker(Marker::Braille) // Wichtig für hohe Auflösung!
        .x_bounds([-10.0, 110.0])
        .y_bounds([-10.0, 110.0])
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
    if let Some(project) = &app.project {
        let points: Vec<(f64, f64)> = project
            .mesh
            .nodes
            .iter()
            .map(|node| project_iso(node.1.x, node.1.y, node.1.z))
            .collect();

        ctx.draw(&Points::new(&points, Color::White));
    }
}

#[allow(clippy::many_single_char_names)]
pub fn project_iso(x: f64, y: f64, z: f64) -> (f64, f64) {
    let cos30 = 0.866;
    let sin30 = 0.5;

    let u = (x - z) * cos30;
    let v = y + (x + z) * sin30;

    (u, v)
}
