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
The project is in its infancy. I am in the progress of implementing the core structure and internal api.
If you still want to try it out I recommend running one of the C3D4 example files in `test_inputs`.
The performance is not far behind what you'd expect from even the commercial solvers.
