use ferrix::solver::ids::{ElementId, NodeId};
use ferrix::solver::mesh_lib::elements::c3d6::C3D6;
use ferrix::solver::mesh_lib::elements::element::FiniteElement;
use nalgebra::DMatrix;

#[test]
fn test_degenerate_wedge_jacobian() {
    // Define a wedge where nodes 2 and 5 are the same (collapsed edge)
    // Bottom triangle: (0,0,0), (1,0,0), (0,1,0)
    // Top triangle:    (0,0,1), (1,0,1), (0,1,1)
    // Collapse 2 and 5 to (1,0,0.5)
    let node_coords = vec![
        [0.0, 0.0, 0.0], // 1
        [1.0, 0.0, 0.5], // 2 (collapsed with 5)
        [0.0, 1.0, 0.0], // 3
        [0.0, 0.0, 1.0], // 4
        [1.0, 0.0, 0.5], // 5 (collapsed with 2)
        [0.0, 1.0, 1.0], // 6
    ];

    let mut coords_mat = DMatrix::<f64>::zeros(3, 6);
    for (i, coord) in node_coords.iter().enumerate() {
        coords_mat[(0, i)] = coord[0];
        coords_mat[(1, i)] = coord[1];
        coords_mat[(2, i)] = coord[2];
    }

    let elem = C3D6 {
        id: ElementId(1),
        nodes: [
            NodeId(1),
            NodeId(2),
            NodeId(3),
            NodeId(4),
            NodeId(5),
            NodeId(6),
        ],
    };

    // Check Jacobian at integration points
    for gp in elem.integration_points() {
        let (_, dn_local) = elem.shape_functions(gp.coords[0], gp.coords[1], gp.coords[2]);
        let jacobian = &dn_local * coords_mat.transpose();
        let det = jacobian.determinant();
        println!("GP at {:?}: det(J) = {}", gp.coords, det);
        assert!(
            det.abs() > 1e-10,
            "Jacobian should not be singular at Gauss point"
        );
    }

    // Check Jacobian at nodes
    let node_local_coords = elem.node_local_coords();
    let centroid = [1.0 / 3.0, 1.0 / 3.0, 0.0];
    for (i, local_pos) in node_local_coords.iter().enumerate() {
        let (_, dn_local) = elem.shape_functions(local_pos[0], local_pos[1], local_pos[2]);
        let jacobian = &dn_local * coords_mat.transpose();
        let det = jacobian.determinant();
        println!("Node {} at {:?}: det(J) = {}", i + 1, local_pos, det);

        if det.abs() < 1e-10 {
            // Try slightly shifted towards centroid
            let eps = 1e-3;
            let shifted = [
                local_pos[0] * (1.0 - eps) + centroid[0] * eps,
                local_pos[1] * (1.0 - eps) + centroid[1] * eps,
                local_pos[2] * (1.0 - eps) + centroid[2] * eps,
            ];
            let (_, dn_local_shifted) = elem.shape_functions(shifted[0], shifted[1], shifted[2]);
            let j_shifted = &dn_local_shifted * coords_mat.transpose();
            let det_shifted = j_shifted.determinant();
            println!(
                "  Shifted Node {} at {:?}: det(J) = {}",
                i + 1,
                shifted,
                det_shifted
            );
            assert!(det_shifted.abs() > 1e-10);
        }
    }
}
