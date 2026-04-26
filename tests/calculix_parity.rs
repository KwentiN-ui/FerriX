use ferrix::solver::{
    ids::NodeId, io::writer::ResultWriter, project::Project, results::IncResult,
    state::SolutionState, time::SolverTime,
};

struct MockWriter;
impl ResultWriter for MockWriter {
    fn init(&self) -> Result<(), Box<dyn std::error::Error>> {
        Ok(())
    }
    fn write_increment(
        &self,
        _res: &IncResult,
        _timer: &SolverTime,
    ) -> Result<(), Box<dyn std::error::Error>> {
        Ok(())
    }
    fn finish(&self) -> Result<(), Box<dyn std::error::Error>> {
        Ok(())
    }
}

#[test]
fn test_calculix_parity_c3d4() {
    let project = Project::from_jobname("tests/data/C3D4_staticstep", None)
        .expect("Failed to parse .inp file");

    let num_dofs = project.mesh.nodes.len() * 3;
    let mut solution_state = SolutionState::new(num_dofs);
    let mut simulation_time = SolverTime::new();
    let writer = MockWriter;

    for (i, step) in project.steps.iter().enumerate() {
        let step_id = i + 1;
        step.solve(
            step_id,
            &project,
            &mut solution_state,
            &writer,
            &mut simulation_time,
        )
        .expect("Step failed");
    }

    // Node 17 Reference from CalculiX .dat file:
    // 17  9.356107E-06 -1.313901E-04  1.719218E-03
    let node_id = NodeId(17);
    let idx = project
        .mesh
        .get_index_for_node_id(node_id)
        .expect("Node 17 not found");

    let ux = solution_state.displacements[idx * 3];
    let uy = solution_state.displacements[idx * 3 + 1];
    let uz = solution_state.displacements[idx * 3 + 2];

    println!("Ferrix results for Node 17: ux={ux:.6e}, uy={uy:.6e}, uz={uz:.6e}");

    let ref_ux = 9.356_107E-06;
    let ref_uy = -1.313_901E-04;
    let ref_uz = 1.719_218E-03;

    // We expect some difference due to solver precision and penalty method
    let tolerance = 1e-6;

    assert!(
        (ux - ref_ux).abs() < tolerance,
        "ux mismatch: ferrix={ux}, ref={ref_ux}"
    );
    assert!(
        (uy - ref_uy).abs() < tolerance,
        "uy mismatch: ferrix={uy}, ref={ref_uy}"
    );
    assert!(
        (uz - ref_uz).abs() < tolerance,
        "uz mismatch: ferrix={uz}, ref={ref_uz}"
    );
}

#[test]
fn test_calculix_parity_2step() {
    let project = Project::from_jobname("tests/data/C3D4_2staticstep", None)
        .expect("Failed to parse .inp file");

    let num_dofs = project.mesh.nodes.len() * 3;
    let mut solution_state = SolutionState::new(num_dofs);
    let mut simulation_time = SolverTime::new();
    let writer = MockWriter;

    // Run Step 1
    project.steps[0]
        .solve(
            1,
            &project,
            &mut solution_state,
            &writer,
            &mut simulation_time,
        )
        .expect("Step 1 failed");

    // Verify Node 17 Step 1
    // 17  9.356107E-06 -1.313901E-04  1.719218E-03
    let node_id = NodeId(17);
    let idx = project
        .mesh
        .get_index_for_node_id(node_id)
        .expect("Node 17 not found");

    {
        let ux = solution_state.displacements[idx * 3];
        let uy = solution_state.displacements[idx * 3 + 1];
        let uz = solution_state.displacements[idx * 3 + 2];

        let ref_ux = 9.356_107E-06;
        let ref_uy = -1.313_901E-04;
        let ref_uz = 1.719_218E-03;
        let tolerance = 1e-6;

        assert!(
            (ux - ref_ux).abs() < tolerance,
            "Step 1 ux mismatch: ferrix={ux}, ref={ref_ux}"
        );
        assert!(
            (uy - ref_uy).abs() < tolerance,
            "Step 1 uy mismatch: ferrix={uy}, ref={ref_uy}"
        );
        assert!(
            (uz - ref_uz).abs() < tolerance,
            "Step 1 uz mismatch: ferrix={uz}, ref={ref_uz}"
        );
    }

    // Run Step 2
    project.steps[1]
        .solve(
            2,
            &project,
            &mut solution_state,
            &writer,
            &mut simulation_time,
        )
        .expect("Step 2 failed");

    // Verify Node 17 Step 2
    // 17  2.603480E-04 -1.377070E-04  1.001299E-03
    {
        let ux = solution_state.displacements[idx * 3];
        let uy = solution_state.displacements[idx * 3 + 1];
        let uz = solution_state.displacements[idx * 3 + 2];

        let ref_ux = 2.603_480E-04;
        let ref_uy = -1.377_070E-04;
        let ref_uz = 1.001_299E-03;
        let tolerance = 1e-6;

        assert!(
            (ux - ref_ux).abs() < tolerance,
            "Step 2 ux mismatch: ferrix={ux}, ref={ref_ux}"
        );
        assert!(
            (uy - ref_uy).abs() < tolerance,
            "Step 2 uy mismatch: ferrix={uy}, ref={ref_uy}"
        );
        assert!(
            (uz - ref_uz).abs() < tolerance,
            "Step 2 uz mismatch: ferrix={uz}, ref={ref_uz}"
        );
    }
}
