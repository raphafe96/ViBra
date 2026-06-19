# ViBra - Vibrational Structure Program

A Fortran program for Vibrational Self-Consistent Field (VSCF) and Vibrational Configuration Interaction (VCI) calculations using ORCA output files.

**______________________________________________________________________**

## Table of Contents

- [Compilation](#compilation)
  - [Prerequisites](#prerequisites)
  - [Compilation Command](#compilation-command)
  - [Compiler Flags](#compiler-flags)
- [Pre-compiled Executable](#pre-compiled-executable)
- [Usage](#usage)
  - [Required Input Files](#required-input-files)
  - [Input Parameters](#input-parameters)
  - [Output Files](#output-files)
- [Program Architecture](#program-architecture)
- [License](#license)

**______________________________________________________________________**

## Compilation

### Prerequisites

- **Intel Fortran Compiler (`ifx`)** — **Required** due to Intel MKL/LAPACK dependencies.
- **Intel MKL** — Included with Intel compiler distribution.
- The code has been tested on **Windows**. Linux/macOS compilation may work but is **untested**.

### Compilation Command

```bash
ifx -O3 /MT /libs:static /Qopenmp /Qm64 /heap-arrays /threads /Qmkl:parallel /fpscomp:logicals main.f90 symmetry.f90 jacobi.f90 integrals.f90 one_mode_operation.f90 read_input.f90 read_orca.f90 get_combination.f90 vci.f90 -o vscf_vci && vscf_vci.exe
```

> **⚠️ Important:** The order of `.f90` files shown above may **not** reflect correct module dependency ordering. Multiple compilation attempts may be necessary, or manually adjust the file order based on module dependencies.

### Compiler Flags

| Flag | Description |
|:-----|:------------|
| `-O3` | Aggressive optimization |
| `/MT` | Static multi-threaded runtime |
| `/libs:static` | Static library linking |
| `/Qopenmp` | OpenMP parallelization |
| `/Qm64` | 64-bit compilation |
| `/heap-arrays` | Heap-based array allocation |
| `/threads` | Multi-threaded MKL |
| `/Qmkl:parallel` | Parallel MKL routines |
| `/fpscomp:logicals` | FPS logical compatibility |

**______________________________________________________________________**

## Pre-compiled Executable

A compiled Windows executable is provided in **`exe.rar`**.

1. Extract the archive.
2. Place the executable together with the required Fortran runtime `.dll` files.

**______________________________________________________________________**

## Usage

### Required Input Files

| File | Description |
|:-----|:------------|
| `input_vscf.txt` | **Main input file.** Must be placed in the same directory as the executable and DLLs. **Do not rename.** |
| `*.vpt2` | ORCA `.vpt2` file containing semi-quartic force field constants. Path must match `FILECT` in `input_vscf.txt`. |

### Input Parameters

**All parameters in `input_vscf.txt` must be set.** Example below:

| Variable | Type | Description |
|:---------|:-----|:------------|
| **NMODES** | `int` | Number of vibrational modes (**M**) |
| **NEXPAN** | `int` | HO basis size per mode (**N**exp) |
| **FILECT** | `str` | Path to ORCA `.vpt2` file |
| **CTEMOD** | `str` | Format specification (`orca_vpt2`) |
| **NQUANT** | `int` | Maximum total quanta (**N**q); ≤ 0 **disables VCI** |
| **NSTATE** | `int` | Number of eigenstates to compute (≤ 0 = **all**) |
| **CVGSCF** | `int` | VSCF convergence exponent (10⁻ᶜ cm⁻¹) |
| **THREAD** | `int` | Number of OpenMP threads |
| **PGROUP** | `str` | Point group (`C1`, `Cs`, `Ci`, `C2`, `C2h`, `C2v`, `D2`, `D2h`) |
| **PROJCT** | `real` | Projection cutoff for symmetry detection (Å) |
| **MAXSCI** | `int` | **N**sel for Selected CI; `0` = full VCI |

### Output Files

| File | Description |
|:-----|:------------|
| **`vscf.out`** | Main output: VSCF/VCI coefficients, symmetry information (if enabled), Selected CI states (if enabled). |
| **`intensities.txt`** | Transition intensities from ORCA `.vpt2` dipole moment derivatives. **1st order** (VSCF/harmonic), **2nd order** (VCI and variants). |
| **`normal_mode.txt`** | Intensities and atomic displacements for each normal mode, including translations and rotations. |

**______________________________________________________________________**

## Program Architecture

```
main.f90                  → Entry point and workflow orchestration
│
├── read_input.f90        → Keyword-based input parser with validation
├── read_orca.f90         → ORCA .vpt2 parser, Hessian diagonalization,
│                            dipole transformation, normal_mode.txt writer
├── symmetry.f90          → Point group setup (character tables, symmetry
│                            operations, direct product tables), normal mode
│                            irrep assignment, block-diagonal VCI dispatch
├── get_combination.f90   → Configuration enumeration (recursive stars-and-bars),
│                            combinatorial counting, unique-element decomposition,
│                            degeneracy factors
├── integrals.f90         → Analytical HO matrix elements (p = 0,…,5)
├── jacobi.f90            → LAPACK wrappers (DSYEVD, DSYEVR), workspace query,
│                            eigenvector normalization
├── one_mode_operation.f90 → VSCF mean-field calculation (OpenMP-parallelized)
└── vci.f90               → Shared compute_H_element kernel, full VCI,
                             Selected VCI with EN-PT2, Symmetry-Adapted VCI,
                             dipole calculations (OpenMP-parallelized)
```

### Module Descriptions

| Module | Description |
|:-------|:------------|
| **`get_combination`** | Configuration enumeration (recursive stars-and-bars), combinatorial counting, unique-element decomposition, degeneracy factors. |
| **`jacobi`** | LAPACK wrappers (DSYEVD, DSYEVR), automatic workspace query, eigenvector normalization. *Named after original Jacobi implementation.* |
| **`integrals`** | Analytical harmonic oscillator matrix elements (powers p = 0,…,5) (p=5 sets kinetic energy). |
| **`read_input`** | Keyword-based input parser with validation. |
| **`read_orca`** | ORCA `.vpt2` parser, Hessian diagonalization, dipole transformation, `normal_mode.txt` writer. |
| **`one_mode_operation`** | VSCF mean-field calculation (OpenMP-parallelized). |
| **`vci`** | Shared `compute_H_element` kernel, full VCI, Selected VCI with EN-PT2, Symmetry-Adapted VCI, dipole calculations (OpenMP-parallelized). |
| **`symmetry`** | Point group setup (character tables, symmetry operations, direct product tables), normal mode irrep assignment via principal axis frame projection, block-diagonal VCI dispatch. |

**______________________________________________________________________**
