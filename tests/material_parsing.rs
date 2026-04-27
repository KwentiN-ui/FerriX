use ferrix::solver::ids::NodeId;
use ferrix::solver::inp::InpFile;
use ferrix::solver::parser::Parser;

#[test]
#[allow(clippy::float_cmp)]
fn test_temperature_dependent_material_parsing() {
    let inp_content = "
*MATERIAL, NAME=TestMat
*DENSITY
1000.0, 0.0
900.0, 500.0
*ELASTIC
210000.0, 0.3, 0.0
150000.0, 0.35, 500.0
";
    let inp = InpFile::new(inp_content);
    let project = Parser::new(&inp).parse().unwrap();

    assert_eq!(project.materials.len(), 1);
    let mat = &project.materials[0];
    assert_eq!(mat.name(), "TESTMAT");

    // Test density interpolation
    assert_eq!(mat.density(0.0), Some(1000.0));
    assert_eq!(mat.density(250.0), Some(950.0));
    assert_eq!(mat.density(500.0), Some(900.0));
    assert_eq!(mat.density(600.0), Some(900.0)); // Clamping

    // Test elastic interpolation
    assert_eq!(mat.youngs_modulus(0.0), Some(210_000.0));
    assert_eq!(mat.poisson_ratio(0.0), Some(0.3));

    assert_eq!(mat.youngs_modulus(250.0), Some(180_000.0));
    let nu_interp = mat.poisson_ratio(250.0).unwrap();
    assert!((nu_interp - 0.325).abs() < 1e-12);
}

#[test]
#[allow(clippy::float_cmp)]
fn test_thermal_expansion_parsing() {
    let inp_content = "
*MATERIAL, NAME=ThermalMat
*EXPANSION, ZERO=20.0
1.2E-5, 0.0
1.5E-5, 500.0
*INITIAL CONDITIONS, TYPE=TEMPERATURE
1, 100.0
";
    let inp = InpFile::new(inp_content);
    let project = Parser::new(&inp).parse().unwrap();

    let mat = &project.materials[0];
    assert_eq!(mat.reference_temperature(), 20.0);
    assert_eq!(mat.thermal_expansion(0.0), Some(1.2E-5));
    assert_eq!(mat.thermal_expansion(250.0), Some(1.35E-5));
    assert_eq!(mat.thermal_expansion(500.0), Some(1.5E-5));

    assert_eq!(project.initial_temperatures.get(&NodeId(1)), Some(&100.0));
}
