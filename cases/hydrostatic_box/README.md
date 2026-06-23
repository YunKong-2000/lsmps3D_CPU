# Hydrostatic Box

This case uses a standalone preprocessor to generate particle input files for a hydrostatic water tank.

Geometry:

- Tank length: 1 m
- Tank height: 0.75 m
- Tank width: 1 m
- Water height: 0.5 m
- Particle spacing: 0.05 m
- Wall representation: single-layer wall particles over the full configured tank walls

Build and generate input particles:

```bash
cmake --build build
./build/hydrostatic_box_preprocess
```

Generated files:

- `cases/hydrostatic_box/fluid_particles.dat`
- `cases/hydrostatic_box/wall_particles.dat`
- `output/hydrostatic_box/preprocess_particles.vtk`

Current generated particle counts:

- Fluid particles: 4851
- Wall particles: 1937
- Total particles: 6788

Run the solver from generated particle files:

```bash
./build/lsmps3d cases/hydrostatic_box/config.ini
```

The configured short validation run outputs:

- `output/hydrostatic_box/hydrostatic_initial.vtk`
- `output/hydrostatic_box/hydrostatic_00000_step_1.vtk`

The time-step VTK output keeps only core fields and selected diagnostics. Detailed debugging fields such as free-surface area ratios, PPE RHS, and body/viscous acceleration vectors are intentionally omitted from this case output.

When a VTK file is written, the solver prints one timing line to the terminal:

```text
[timing] step=... neighbor=... free_surface=... lsmps=... provisional=... ppe=... correction=... vtk=... total=...
```

Particle file format:

- Fluid rows: `x y z vx vy vz pressure density`
- Wall rows: `x y z vx vy vz pressure density nx ny nz`
