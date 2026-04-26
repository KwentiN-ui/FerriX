# CalculiX Solver Structure and Solution Process Report

_This report is **AI generated**, and highlights the core calculix code-structure._

This report outlines the architecture and execution flow of the CalculiX Solver (CCX), with a focus on Mechanical and Thermal analysis and the handling of subsequent steps.

## 1. Overall Codebase Structure

CalculiX is a hybrid C and Fortran program. The high-level orchestration and memory management are primarily handled in C, while the heavy numerical computations (element assembly, material models, results calculation) are implemented in Fortran.

### Key Components:
- **Entry Point:** `src/CalculiX.c` contains the `main` function.
- **Input Parsing:** `src/readinput.c` (initial read) and `src/calinput.f` (step-by-step reading).
- **Memory Management:** `src/allocation.f` handles the large-scale array allocations.
- **Solvers:**
  - `src/nonlingeo.c`: The primary solver for non-linear geometry and coupled mechanical/thermal problems.
  - `src/linstatic.c`: Linear static analysis.
  - `src/arpack.c` / `src/arpackcs.c`: Frequency analysis using ARPACK.
- **Assembly:** `src/mafillsmmain.c` (C wrapper for multithreading) and `src/mafillsm.f` (Fortran assembly logic).
- **Element Calculations:** `src/e_c3d.f` (and variants) for 3D elements.
- **Material Models:** `src/mechmodel.f` and `src/materialdata_me.f`.
- **Results & Internal Forces:** `src/results.c` (C wrapper) and `src/resultsmech.f` (Fortran logic).

---

## 2. Call Hierarchy and Solution Process

The execution flow follows a **Research -> Strategy -> Execution** pattern internally:

1.  **Initialization:** `main` (CalculiX.c) reads the job name, opens files, and calls `readinput` to parse the model definition (nodes, elements, materials).
2.  **Step Loop:** A `while(istat>=0)` loop in `main` processes each `*STEP` defined in the input deck.
3.  **Step Input:** `calinput.f` is called to read the current step's parameters (analysis type, loads, boundary conditions).
4.  **Analysis Execution:** Based on the `nmethod` variable (determined in `calinput`), the appropriate solver is called. For mechanical/thermal problems, this is usually `nonlingeo`.

### The Non-Linear Solution Loop (`nonlingeo.c`):
`nonlingeo` implements a Newton-Raphson iteration scheme within an increment loop:

- **Increment Loop:** Divides the step time into smaller increments.
- **Iteration Loop (Newton-Raphson):**
  - **Assembly:** Calls `mafillsmmain` -> `mafillsm.f`. This calculates the Tangent Stiffness Matrix ($K$).
    - `mafillsm.f` loops over elements and calls `e_c3d.f`.
    - `e_c3d.f` performs numerical integration (Gauss) and calls `materialdata_me.f` for constitutive properties.
  - **Residual Calculation:** Calls `calcresidual.c`. It computes the out-of-balance forces $R = F_{ext} - F_{int}$.
  - **Internal Forces:** `results` -> `resultsmech.f` computes $F_{int}$ by integrating stresses over the elements.
    - `resultsmech.f` calls `mechmodel.f` to get stresses from strains/state variables.
  - **Linear Solve:** Calls an external solver (e.g., SPOOLES, PARDISO) to solve $K \Delta u = R$.
  - **Update:** Updates the displacement field $u = u + \Delta u$.
  - **Convergence Check:** Checks if $R$ and $\Delta u$ are within tolerances.

---

## 3. Mechanical and Thermal Analysis

### Mechanical Analysis:
- Focuses on solving the equilibrium equations.
- Handles non-linearities:
  - **Geometric (NLGEOM):** Updated Lagrangian formulation. Strains are calculated in `calctotstrain.f`.
  - **Material:** Plasticity, creep, hyperelasticity. Managed in `mechmodel.f`.
  - **Contact:** Handled via penalty or mortar methods in `contact.f`.

### Thermal Analysis:
- Solves the heat conduction equation.
- Can be **Uncoupled** (pure heat transfer) or **Coupled** (thermo-mechanical).
- Managed by `ithermal` variable (1: Mech, 2: Therm, 3: Coupled).
- Element thermal matrices are computed in `e_c3d_th.f`.

---

## 4. Handling Subsequent Steps

CalculiX supports multi-step simulations where the state of the model at the end of one step is the starting point for the next.

1.  **State Preservation:** The `main` loop ensures that variables like displacements (`v`), stresses (`sti`), and state variables (`xstate`) are preserved across steps.
2.  **Load/BC Updates:** `calinput.f` reads the new step's `*BOUNDARY` and `*CLOAD`/`*DLOAD` cards. It can either add to existing loads or replace them (controlled by the `OP` parameter).
3.  **Step-Specific Logic:** Some analysis types (like `*FREQUENCY`) can be performed using the pre-stress state from a previous `*STATIC` step.
4.  **Re-Allocation:** If a step introduces new constraints or contact pairs, `main` may call `allocation` again or perform re-allocation to accommodate the new matrix structure.

---

## 5. Summary for Rust Implementation

To build a drop-in replacement in Rust, the following architectural layers are recommended:
1.  **Input Parser:** A robust parser for Abaqus-style `.inp` files (similar to `readinput` and `calinput`).
2.  **Finite Element Library:** Supporting various integration schemes and element topologies (similar to `e_c3d`).
3.  **Constitutive Model Interface:** A modular system for material models (similar to `mechmodel`).
4.  **Global Assembler:** A multithreaded sparse matrix assembler (similar to `mafillsmmain`).
5.  **Non-Linear Solver:** A Newton-Raphson orchestrator with time-stepping logic (similar to `nonlingeo`).
6.  **Linear Solver Interface:** Integration with sparse solvers like Faer or others available in the Rust ecosystem.

---

## 6. Current Rust (Ferrix) Implementation vs. CalculiX

### Step Handling Analysis
The current Ferrix implementation handles multiple steps through a loop in `main.rs`, passing a mutable `SolutionState` between steps. However, several architectural differences exist compared to CalculiX:

1.  **Load & Boundary Condition Management:**
    *   **Ferrix:** Loads and BCs are parsed into global vectors (`project.loads`, `project.bcs`). The solver applies ALL loads in these vectors during every step. There is no concept of step-specific scoping or the `OP=NEW/MOD` parameter found in CCX.
    *   **CalculiX:** `calinput.f` reads loads/BCs step-by-step. By default (`OP=MOD`), it modifies existing ones; with `OP=NEW`, it clears previous ones. This ensures only relevant loads are active.

2.  **Displacement Accumulation (Linear Static):**
    *   **Ferrix:** The `StaticStep::solve` method adds the displacement result of the current step (`delta_u`) to the existing `solution_state.displacements`. Since `delta_u` is calculated using the total force (including previous steps' loads), this results in double-counting displacements.
    *   **CalculiX:** For linear static analysis (`linstatic.f`), CCX solves for the total state. For non-linear analysis (`nonlingeo.c`), it solves for increments using residuals ( = F_{ext} - F_{int}$), which naturally handles accumulation.

3.  **Amplitude & Time Ramping:**
    *   **Ferrix:** The default amplitude (`None`) ramps a load from 0 to 1 relative to the *current step's* time. If a load from Step 1 persists in Step 2, Ferrix will incorrectly ramp it from 0 again instead of holding its final value.
    *   **CalculiX:** Loads without an explicit `*AMPLITUDE` are ramped in the step they are introduced and then held constant in subsequent steps.

4.  **Parsing Scoping:**
    *   **Ferrix:** The `Parser` does not currently distinguish between global data and step-specific data. Keywords like `*CLOAD` and `*BOUNDARY` are processed identically regardless of whether they appear inside a `*STEP` block.
    *   **CalculiX:** The parser is context-aware. Data before the first `*STEP` is global (mostly model definition), while data within `*STEP` blocks is local to that analysis phase.

### Recommendations for CCX-Parity:
*   **Encapsulate Step Data:** Move `loads` and `bcs` from the global `Project` struct into the `Step` enum/structs, or implement a stateful manager that handles `OP=NEW/MOD`.
*   **Incremental Force Formulation:** Transition the solver to use an incremental approach ( \Delta u = \Delta F$ or  \Delta u = F_{ext} - F_{int}$) to ensure correct state accumulation across steps.
*   **Improve Amplitude Logic:** Update `Amplitude` to handle "Step-Relative" vs. "Total-Simulation" time correctly, ensuring loads persist across steps as expected.

---

## 7. Refactoring Progress (Rust Implementation)

The following improvements have been implemented in the Rust version to align more closely with CalculiX's behavior:

1.  **Step-Specific Data Scoping:**
    *   The `Parser` now distinguishes between global (initial) loads/BCs and step-specific ones.
    *   `StaticStep` now holds its own collections of `loads` and `bcs`.
2.  **Correct Displacement Accumulation:**
    *   The solver now calculates the **total** displacement for each increment/step using the total active force.
    *   `SolutionState` is updated with the total state at the end of each step, preventing the previous double-counting error.
3.  **Improved Amplitude and Load Persistence:**
    *   Introduced `origin_step` tracking for all loads and BCs.
    *   Updated `Amplitude` logic: Loads without explicit amplitude definitions now stay at their final value (1.0) in subsequent steps instead of re-ramping from zero.
4.  **Parser Context Awareness:**
    *   The `Parser` now uses a `step_counter` and tracks the `current_step` context, correctly assigning data to the active step or the global project state.

**Current Status:** The core architecture now supports basic multi-step analysis with correct state accumulation and load persistence, matching the fundamental logic of CalculiX's linear static steps.

---

## 8. Nonlinear Analysis Handling

Currently, the Rust implementation does **not** support nonlinear analysis. 

### Status:
1.  **Linear Formulation:** The `StaticStep::solve` and `StaticStep::next_increment` methods use a linear static formulation. They assemble the stiffness matrix once per increment based on the initial geometry and solve for the total displacement.
2.  **Newton-Raphson Loop:** There is no iterative loop to minimize residuals ( = F_{ext} - F_{int}$) within an increment.
3.  **Constitutive Models:** `Material` is currently restricted to linear elasticity. Nonlinear material models (plasticity, hyperelasticity) are not yet implemented.
4.  **Geometric Nonlinearity:** Element stiffness calculations (`compute_stiffness`) use the initial configuration. An updated Lagrangian formulation for `NLGEOM` support is missing.

### Recommended Path for Nonlinear Support:
*   **Newton-Raphson Orchestrator:** Implement an iteration loop in `StaticStep` that checks for force and displacement convergence.
*   **Residual Calculation:** Implement a method to calculate the internal force vector {int}$ by integrating stresses over element volumes.
*   **Geometric Nonlinearity:** Add support for tracking the deformed configuration and updating the Jacobian/B-matrix accordingly.
*   **Nonlinear Material Interface:** Extend the `Material` trait to handle state-dependent constitutive updates (e.g., return tangent stiffness and updated stress).

# Roadmap for contacts
1.  **Refactor `StaticStep` to Non-Linear:** 
    *   Implement an inner Newton-Raphson iteration loop.
    *   Implement the calculation of the global Internal Force vector ($F_{int}$).
    *   Add residual convergence checking ($||F_{ext} - F_{int}|| < tol$).
2.  **Geometric Nonlinearity (NLGEOM):** 
    *   Implement updated coordinates ( deformed mesh ).
    *   Calculate stiffness based on the deformed state.
3.  **Basic Constraint Equations (MPCs):** 
    *   Implement Multi-Point Constraints (`*EQUATION`). This lays the mathematical groundwork for tying nodes together, which is a simpler precursor to tying them together conditionally (contact).
4.  **Surface Projections:** 
    *   Build the math to project a point onto a 3D triangle/quadrilateral.
5.  **Finally, Contact Mechanics:** 
    *   Parse `*CONTACT PAIR`.
    *   Implement a basic Node-to-Surface Penalty contact algorithm using the Newton-Raphson loop built in Step 1.
