# Next Steps: Solver Evolution

With the material trait refactoring, temperature-dependent properties, thermal expansion, and SDV infrastructure implemented, the following steps are proposed to enhance the solver's physical capabilities and performance.

## 1. Non-Linear Material Models (Plasticity)
Now that the SDV (Solution-Dependent State Variable) infrastructure is in place, we can implement path-dependent materials.
- **J2 Plasticity:** Implement a `PlasticMaterial` struct using von Mises yield criteria and isotropic hardening.
- **Radial Return Mapping:** Implement the return mapping algorithm within the `update_state` method to calculate updated stress and tangent stiffness.
- **Verification:** Verify against CalculiX's `*PLASTIC` behavior.

## 2. Prescribed Temperature Evolution
We have a global temperature field, but it is currently static (initial conditions only).
- **`*TEMPERATURE` Support:** Update the `Parser` to handle temperature changes within a step.
- **Time-Interpolation:** Update `SolutionState` to interpolate nodal temperatures over time during an increment, allowing for realistic thermal stress analysis during transient steps.

## 3. Distributed Loads (`*DLOAD`)
Currently, only concentrated nodal loads are supported.
- **Pressure Loading:** Implement the `*DLOAD` card for applying uniform pressure to element faces.
- **Surface Integration:** Update `Element` logic to perform surface integration for equivalent nodal forces.

## 4. Higher-Order Elements (`C3D10`)
Linear tetrahedrons (`C3D4`) are prone to shear locking and require very fine meshes.
- **Second-Order Tetra:** Implement the `C3D10` element.
- **Shape Functions:** Add quadratic shape functions and 4-point or 5-point Gauss integration schemes.

## 5. Parallel Assembly
The current `Assembler` iterates through elements sequentially.
- **Rayon Integration:** Use `rayon` to parallelize the element loop in `Assembler::assemble` and `assemble_internal_force`.
- **Thread-Safe Accumulation:** Use a thread-safe sparse matrix builder or a coloring/reduction strategy to assemble the global system.

## 6. Orthotropic & Anisotropic Support
- **Material Symmetries:** Implement `OrthotropicMaterial` (9 constants) and `AnisotropicMaterial` (21 constants).
- **Parser Intelligence:** Enhance the `Parser` to automatically detect the material symmetry based on the number of constants provided in the `*ELASTIC` card.
