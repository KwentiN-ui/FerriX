# FerriX
This project aims to recreate the [CalculiX Solver](https://www.calculix.de/) by Guido Dhondt, meant for educational purposes.
It aims to work with the same input file structure, so it can be used as a drop-in replacement for a small subset of supported problems.
Being written fully in Rust, it allows for more intuitive abstractions than Fortran, while preserving good performance.

The codebase is meant to be as modular as possible, allowing everyone to fork the project and try out their own elements/solvers/steptypes.
I intent to fully leverage the Rust Typesystem to make the code as straightforward as possible. Given the scope of the project I am currently
using AI (mostly Gemini) to iterate fast on the prototype. Nonetheless there is no vibecoding involved.

<img width="517" height="537" alt="impeller" src="https://github.com/user-attachments/assets/e7f470ac-64c7-486f-bbfb-4164abba4c90" />

_Impeller Model solved with FerriX and postprocessed in Paraview_

## Installation

See the releases tab to download prebuilt binaries or use
```
cargo install --git https://github.com/KwentiN-ui/ferrix.git
```
to compile the newest version yourself. This will also make the tool available as `ferrix` in your terminal.

## Features

FerriX aims for full numerical parity with CalculiX for all supported features. Supported keywords include:

### Analysis Steps
- `*STATIC` (Linear and `NLGEOM`)

### Element Types
- `C3D4` (4-node linear tetrahedron)
- `C3D20` (20-node quadratic brick)

### Material Definition
- `*MATERIAL`
- `*ELASTIC` (Temperature-dependent)
- `*DENSITY` (Temperature-dependent)
- `*EXPANSION` (Temperature-dependent, with `ZERO`)
- `*DEPVAR` (Solution-dependent state variables)

### Boundary Conditions & Loading
- `*BOUNDARY`
- `*CLOAD`
- `*INITIAL CONDITIONS, TYPE=TEMPERATURE`
- `*AMPLITUDE`

### Mesh & Sets
- `*NODE`
- `*ELEMENT`
- `*NSET`
- `*ELSET`
- `*SOLID SECTION`

### Output Control
- `*NODE FILE` (U, RF)
- `*EL FILE` (S, E)
- `*NODE PRINT`
