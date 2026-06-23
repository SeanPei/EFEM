# EFEM

An efficient nonlinear finite element solver for isotropic hyperelastic materials — L-BFGS with an approximate Jacobian preconditioner based on the total Lagrangian formulation.

Reference: Pei, X., Du, J., & Chen, G. (2023). *Accelerated nonlinear finite element method for analysis of isotropic hyperelastic materials nonlinear deformations*. Applied Mathematical Modelling, 120, 513–534.

## Features

- Total Lagrangian FEM for large-strain hyperelastic analysis
- Approximate Jacobian matrix that avoids repeated assembly and factorization
- L-BFGS solver with the approximate Jacobian as the initial Hessian estimate
- 3 hyperelastic material models: Neo-Hookean, Mooney-Rivlin, Yeoh
- 6 element types with full and reduced integration: C3D4, C3D8, C3D10, C3D20
- External loads: fixed boundary, nodal force, surface pressure

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| [VTK](https://vtk.org/) | ≥ 8.0 | Mesh I/O and visualization |
| [Boost](https://www.boost.org/) | ≥ 1.65 | Filesystem and string utilities |
| [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) | 2022.1.0 | LAPACK/BLAS for Eigen |
| [Eigen](https://eigen.tuxfamily.org/) | 3.3.7 (bundled) | Linear algebra |
| [spdlog](https://github.com/gabime/spdlog) | 1.x (bundled) | Logging |
| [nlohmann/json](https://github.com/nlohmann/json) | 3.x (bundled) | JSON parsing |
| [CLI11](https://github.com/CLIUtils/CLI11) | 2.x (bundled) | Command-line parsing |

Bundled libraries are in `external/`. VTK, Boost, and MKL must be installed separately.

## Build

```bash
# Configure (adjust paths to your MKL/VTK/Boost installations)
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build --config Release
```

## Usage

```bash
EFEM.exe -f <config.json>
```

Paths in the config file use `/` as the separator on all platforms.

## Configuration

The input is a JSON file with the following sections:

```json
{
  "analysis": {
    "max_optim_iters": 0,
    "max_line_search": 20
  },
  "model": [
    {
      "type": "soft",
      "mesh_file": "mesh.txt"
    }
  ],
  "constraint": [
    { "id": 1, "type": "fix",     "constraint_file": "fix.txt" },
    { "id": 2, "type": "pressure", "value": 10000, "constraint_file": "pressure.txt" }
  ],
  "material": [
    {
      "id": 1,
      "type": "neo_hookean",
      "density": 1050.0,
      "C10": 114000,
      "D1": 1e-6
    }
  ],
  "output": {
    "path": "results"
  }
}
```

### `analysis`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `max_optim_iters` | int | `0` | Max L-BFGS iterations (`0` = unlimited) |
| `max_line_search` | int | `20` | Max backtracking line-search trials |
| `dx_criterion` | double | `1e-8` | Displacement convergence tolerance |
| `residual_criterion` | double | `1e-6` | Residual convergence tolerance |

### `model`

| Key | Type | Description |
|-----|------|-------------|
| `type` | string | `"soft"` for deformable bodies |
| `mesh_file` | string | Path to mesh file (Gmsh format, `.msh` extension auto-detected) |

Mesh file format: one line per node or element. Nodes: `x y z`. Elements: `node_i node_j node_k ...` (indices are 1-based).

### `constraint`

| Type | Key | Description |
|------|-----|-------------|
| `fix` | `constraint_file` | File listing node indices to fix (comma or space separated) |
| `force` | `value: [fx, fy, fz]`, `constraint_file` | Uniform nodal force on listed nodes |
| `pressure` | `value: p` or `value: [px, py, pz]`, `constraint_file` | Scalar pressure or directed pressure on surface |

### `material`

| Type | Parameters |
|------|------------|
| `"neo_hookean"` | `C10`, `D1`, `density` |
| `"mo_riv"` | `C10`, `C01`, `D1`, `density` |
| `"yeoh"` | `C10`, `C20`, `C30`, `D1`, `density` |

## Examples

The `sample/` directory contains the three test cases from the paper:

| Example | Element variants | Description |
|---------|-----------------|-------------|
| `cap/` | C3D4, C3D20 | Rubber sealing cap under pressure |
| `beam/` | C3D4, C3D20 | Cantilever beam under surface traction |
| `finger/` | C3D4, C3D10 | Pneumatic-driven soft finger |

Run with:

```bash
EFEM.exe -f sample/cap/cap.json
EFEM.exe -f sample/beam/beam.json
EFEM.exe -f sample/finger/finger.json
```

Edit the `"mesh_file"` and `"constraint_file"` paths in each JSON to switch between element types.

## Output

The solver writes deformed meshes as VTK files (`undeformed*.vtk`, `deformed*.vtk`) that can be opened in ParaView.

## License

MIT License. See [LICENSE.txt](LICENSE.txt).
