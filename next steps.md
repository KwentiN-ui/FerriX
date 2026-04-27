# Next Steps: Material System & Thermal Coupling

Following the successful refactoring of the `Material` system into a trait-based model with temperature-dependent Look-Up Tables (LUTs), the following steps are proposed to further align FerriX with CalculiX capabilities and Abaqus-style APIs.

## 1. Thermal Expansion & Reference Temperatures
To support thermal loading, we need to implement thermal strain calculations.
- **Trait Extension:** Add `thermal_expansion(&self, temp) -> Option<f64>` to the `Material` trait.
- **Reference State:** Add `reference_temperature` to the `Project` or `BaseMaterial` (Abaqus `*EXPANSION, ZERO=...`).
- **Solver Integration:** Update the `Assembler` to calculate the thermal strain vector $\epsilon_{th} = \alpha(T)(T - T_{ref})$ and add it to the internal force vector or load vector.

## 2. Solution-Dependent State Variables (SDVs)
For non-linear materials (plasticity, creep), the material needs to "remember" its state.
- **State Tracking:** Create a `MaterialState` struct to hold state variables (similar to CCX `STAT` or Abaqus `SDVs`).
- **Trait Evolution:** Update `Material` methods to optionally accept/return updated state variables.
- **Storage:** Implement a mechanism in `Project` or `SolutionState` to store these variables at each integration point.

## 3. Orthotropic & Anisotropic Support
The current `BaseMaterial` assumes isotropy.
- **New Implementations:** Create `OrthotropicMaterial` and `AnisotropicMaterial` structs implementing the `Material` trait.
- **Matrix Logic:** Implement the specialized `build_elastic_d_matrix` for these symmetry classes (9 constants for orthotropic, 21 for fully anisotropic).
- **Parser Extension:** Update `Parser::parse_elastic` to detect the number of constants and instantiate the correct material type.

## 4. User-Defined Materials (Rust "UMAT")
Leverage the trait-based design to allow library users to define materials in external crates.
- **Documentation:** Provide an example of a custom struct implementing the `Material` trait.
- **Plugin System:** Explore the possibility of loading shared libraries (.so/.dll) that implement a C-ABI version of the `Material` trait, allowing users to drop in existing Fortran UMATs.

## 5. Global Temperature Field
Currently, the solver passes a hardcoded `0.0` to material calls.
- **State Field:** Add a `temperatures: Vec<f64>` field to `SolutionState` (nodal temperatures).
- **Interpolation:** Use element shape functions to interpolate nodal temperatures to integration points before calling material methods.
- **Initial Conditions:** Add support for `*INITIAL CONDITIONS, TYPE=TEMPERATURE` in the parser.
