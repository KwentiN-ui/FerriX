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
        .constraints([Constraint::Percentage(70), Constraint::Min(30)])
        .split(f.area());

    // Mesh canvas
    let canvas = Canvas::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Mesh Preview ")
                .title_bottom("Zoom - X/Y"),
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
    if let Some(project) = &app.project {
        // Calculate MVP matrix once per frame
        let mvp = app.camera.build_view_projection_matrix();

        let points: Vec<(f64, f64)> = project
            .mesh
            .nodes
            .iter()
            .filter_map(|node| {
                let p_world = Point3::new(node.1.x, node.1.y, node.1.z);

                // transform_point performs perspective division automatically
                let p_ndc = mvp.transform_point(&p_world);

                // Frustum culling: check if point is within NDC bounds [-1.0, 1.0]
                if p_ndc.x.abs() <= 1.0 && p_ndc.y.abs() <= 1.0 && p_ndc.z.abs() <= 1.0 {
                    Some((p_ndc.x, p_ndc.y))
                } else {
                    None
                }
            })
            .collect();

        ctx.draw(&Points::new(&points, Color::White));
    }
}
