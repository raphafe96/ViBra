<table border="0">
  <tr>
    <td><img src="assets/logo.png" alt="ViBra Logo" width="300"/></td>
    <td><strong>A Selected and Symmetry-Adapted VCI Platform for Anharmonic
Vibrational Spectroscopy and Quantum-Sampled Configuration Spaces.</strong><br><br><em>The logo is inspired on the sunset at Morro Dois Irm&atilde;os, a classic landscape postcard from Rio de Janeiro, resembling two overlapping Gaussian peaks.</em></td>
  </tr>
</table>

> **Citation**: If you use ViBra in your research, please refer to:  
> **(https://doi.org/10.48550/arXiv.2607.22850)**  
>  
> **Theory and Manual**: [📖 Read (English)](manual.pdf) | [📖 Ler (Português)](manual_pt.pdf)
## ✨ Main Features

* Harmonic oscillator, VSCF, VCI, Selected VCI, and Symmetry-Adapted VCI calculations
* Support for cubic and quartic force fields obtained from ORCA VPT2 calculations
* Infrared intensities using first- and second-order dipole derivatives
* OpenMP parallelization
* LAPACK/BLAS diagonalization routines
* Symmetry treatment for the Abelian point groups
* Output files compatible with the ViBra graphical spectrum viewer

## 🔬 Playground

New and experimental features are being developed in the [`playground`](playground/) folder. Check there to see what is new or different from the stable version.

## 📁 Repository Structure

    ViBra/
    │
    ├── source/
    │   ├── main.f90
    │   ├── read_input.f90
    │   ├── read_orca.f90
    │   ├── get_combination.f90
    │   ├── integrals.f90
    │   ├── one_mode_operation.f90
    │   ├── vci.f90
    │   ├── jacobi.f90
    │   ├── symmetry.f90
    │
    ├── playground/
    │   ├── new versions of the code with new features
    │
    ├── examples/
    │   ├── examples.zip
    │
    ├── GUI_and_precompiled_Windows/
    │   ├── Final_dist.zip
    │
    ├── assets/
    │   └── logo.png
    │
    ├── README.md
    └── LICENSE

## ⚠️ Important Note About the GUI and Precompiled Executable

The graphical user interface and the precompiled Fortran executable are provided in a separate folder:

    GUI_and_precompiled_Windows/

These files are intended to work only on Windows.

The precompiled executable was built for Windows, and the GUI was developed for use on Windows systems. Users working on Linux or macOS should compile the Fortran source code themselves and run the program from the command line.

**All files should be placed in the same directory. Running the GUI (vscf_vci_gui.exe) will automatically copy the Fortran executable (vscf_vci.exe) and all dependencies into the user-selected input file directory, and erase them after execution.**

## ⚠️ A Note About ORCA Files

To generate the required `.vpt2` file for ViBra, you need to run two sequential ORCA calculations. Below are examples using **ORCA 6.1**.

**1. Geometry Optimization**

    ! Opt VeryTightSCF wB97X-D4 aug-cc-pvtz

    %scf
       MaxIter 300
    end

    %geom
       MaxIter 300
       Calc_Hess true
       Recalc_Hess 10
       TolE 1e-12
       TolRMSG 1e-8
       TolMaxG 1e-8
       TolRMSD 1e-8
       TolMaxD 1e-8
    end

    * xyz 0 1
      [your geometry here]
    *

**2. VPT2 Calculation**

    ! VeryTightSCF wB97X-D4 aug-cc-pvtz VPT2

    %vpt2
      VPT2 true
      PrintLevel 2
    end

    * xyz 0 1
      [your optimized geometry here]
    *

This will produce a `basename.vpt2` file containing harmonic frequencies, normal modes, cubic and quartic force constants, and first- and second-order dipole derivatives — all required by ViBra.

**Note:** These calculations were performed and tested using **ORCA 6.1**. For more details, see:

> Neese, F. *et al.* ORCA – An Ab Initio, DFT and Semiempirical SCF-MO Package, Version 6.1. Max-Planck-Institut für Kohlenforschung, Mülheim an der Ruhr, 2025. Available at: https://www.faccts.de/orca

## 🔧 Requirements for Compilation

To compile ViBra from source, the following software is required:

* A Fortran compiler with Fortran 90/95 support
* OpenMP support
* LAPACK
* BLAS

The code was developed and tested using Intel Fortran (`ifx`). Other compilers may work, but the compilation flags and linked libraries may need to be adjusted.

Recommended compiler:

    Intel Fortran Compiler (ifx)

Recommended numerical libraries:

    Intel Math Kernel Library (MKL)

> **Note for larger molecules**: When compiling with `ifx`, it may be necessary to use the `-i8` flag to enable 64-bit integers and to set a fixed stack size for the binary. This can help avoid memory addressing issues for very large VCI spaces. Consult your compiler documentation for the appropriate flags.

## 🖥️ Compilation on Windows Using Intel Fortran

An example compilation command using Intel Fortran is:

    ifx /O3 /Qopenmp /threads /Qmkl:parallel /Qm64 /heap-arrays /fpscomp:logicals *.f90 /exe:ViBra.exe

This command enables optimization, OpenMP parallelization, threaded MKL routines, 64-bit compilation, and heap allocation for temporary arrays.

> **⚠️ Important Note — Non-MKL Compilation and File Order**:
>
> **If you are not using Intel MKL**, you must edit `main.f90` before compiling. The Intel MKL dependency appears in two places within that file:
>
> 1. **Module declaration** — remove the line:
>    ```fortran
>    use mkl_service
>    ```
> 2. **Thread setting** — remove or replace the line:
>    ```fortran
>    call mkl_set_num_threads(...)
>    ```
>    with an equivalent call from your LAPACK/BLAS library (e.g., `openblas_set_num_threads` for OpenBLAS). **There are two occurrences of this declaration in main.f90. Replace both!**
>
> All other source files are portable: `jacobi.f90` uses standard LAPACK routines (`dsyevr`, `dsyevd`), and `vci.f90` uses `ddot` from BLAS. You only need to ensure that your compiler links against an available LAPACK/BLAS wrapper (e.g., `-llapack -lblas` for reference implementations, or `-lopenblas` for OpenBLAS).
>
> **Compilation order matters.** Each `.f90` file declares its own modules at the very beginning, and later files depend on modules defined in earlier files. Using a wildcard like `*.f90` may not resolve dependencies correctly and can lead to incorrect results. Compile the files module-by-module, checking the dependencies declared at the top of each `.f90` file to determine the correct order. Alternatively, you can force compile multiple times until all module dependencies are resolved.

If additional modules are included in the source directory, compile them before the files that use them.

## 🐧 Compilation on Linux

A possible compilation command using Intel Fortran on Linux is:

    ifx -O3 -qopenmp -qmkl=parallel -heap-arrays *.f90 -o ViBra

If using `gfortran`, LAPACK, BLAS, and OpenMP must be linked manually. For example:

    gfortran -O3 -fopenmp *.f90 -llapack -lblas -o ViBra

The exact command may vary depending on the installed libraries and operating system.

> **⚠️ Important Note for Non-MKL Users**: The same MKL-dependency and compilation-order warnings described above apply here as well. If you are not using MKL, edit `main.f90` to remove the `use mkl_service` statement and the `call mkl_set_num_threads(...)` line, replacing the latter with an equivalent thread-setting call from your LAPACK/BLAS library. Compile the files module-by-module by checking the dependencies at the top of each `.f90` file, or force compile multiple times until all module dependencies are resolved.

## 📝 Input File

ViBra reads the input file:

    input_vscf.txt

The input file uses keyword-based entries. A typical example is:

    NMODES 12
    NEXPAN 10
    FILECT molecule.vpt2
    CTEMOD orca_vpt2
    NQUANT 4
    NSTATE 20
    CVGSCF 6
    THREAD 8
    PGROUP D2h
    PROJCT 0.01
    MAXSCI 0

## 📋 Input Keywords

| Keyword  |    Type | Description                                                                                                                                                                                                                                                                                            |
| -------- | ------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `NMODES` | Integer | Number of vibrational modes (M).                                                                                                                                                                                                                                                               |
| `NEXPAN` | Integer | HO basis size per mode (Nexp).                                                                                                                                                                                                                                                                  |
| `FILECT` |  String | Path to ORCA `.vpt2` file.                                                                                                                                                                                                                                                                      |
| `CTEMOD` |  String | Format (orca_vpt2).                                                                                                                                                                                                                                                                              |
| `NQUANT` | Integer | Maximum total quanta (Nq); ≤ 0 disables VCI.                                                                                                                                                                                                                                                    |
| `NSTATE` | Integer | Number of eigenstates (≤ 0 = all).                                                                                                                                                                                                                                                               |
| `CVGSCF` | Integer | VSCF convergence exponent (10−C cm−1).                                                                                                                                                                                                                                                          |
| `THREAD` | Integer | Number of OpenMP threads.                                                                                                                                                                                                                                                                       |
| `PGROUP` |  String | Point group (C1, Cs, Ci, C2, C2h, C2v, D2, D2h).                                                                                                                                                                                                                                                       |
| `PROJCT` |    Real | Projection cutoff for symmetry detection (Å).                                                                                                                                                                                                                                                   |
| `MAXSCI` | Integer | Nsel for S-VCI. Default: 100.<br><br>**Usage:** `MAXSCI [N] [mode] [ref]`<br>N > 0, mode = auto: iterative selection with CI reference (ref = s, d, t, q).<br>N > 0, mode = list: iterative selection with list reference, no ref needed.<br>N = 0: full VCI (no S-VCI), no mode and no ref needed.<br>N = 0, mode = list: full VCI using user-provided state list, no ref needed. |

## ▶️ Running ViBra

After compiling the code, place the executable in the same directory as `input_vscf.txt`, or provide the correct path to the input file and ORCA output file.

On Windows:

    ViBra.exe

On Linux or macOS:

    ./ViBra

The program reads `input_vscf.txt` and starts the calculation.

## 📥 ORCA Input Requirement

ViBra requires an ORCA `.vpt2` output file containing the anharmonic vibrational information, including:

* Harmonic frequencies
* Normal modes
* Cubic force constants
* Quartic force constants
* First-order dipole derivatives
* Second-order dipole derivatives, when available

The ORCA calculation must be configured to generate the required VPT2 output.

## 📤 Output Files

ViBra produces the following main output files:

### `vscf.out`

This file contains the complete calculation log, including:

* VSCF convergence information
* VSCF modal coefficients
* VCI configuration list
* Vibrational energies
* VCI eigenvectors
* Leading configuration interaction coefficients
* Vibrational state assignments

### `intensities.txt`

This file contains vibrational transition frequencies and normalized infrared intensities for:

* Harmonic oscillator calculations
* VSCF calculations
* VCI, S-VCI, and SA-VCI calculations

### `normal_mode.txt`

This file contains:

* Equilibrium molecular geometry
* Cartesian normal-mode displacement vectors

This file can be used by the graphical viewer to visualize vibrational normal modes.

## 🧮 Calculation Modes

### Harmonic Oscillator Calculation

Set:

    NQUANT 0

This disables the VCI calculation and performs the harmonic and VSCF-related steps.

### Full VCI Calculation

Set:

    MAXSCI 0

and choose a positive value for `NQUANT`.

Example:

    NQUANT 4
    MAXSCI 0

### Selected VCI Calculation

Set `MAXSCI` to a positive value.

Example:

    NQUANT 6
    MAXSCI 100 auto d

A larger `MAXSCI` value retains more configurations and generally improves agreement with full VCI, at the cost of additional computational time.

### Symmetry-Adapted VCI Calculation

Set the molecular point group using `PGROUP`.

Example:

    PGROUP D2h

Symmetry-adapted VCI reduces computational cost by block-diagonalizing the VCI Hamiltonian according to irreducible representations.

The selected point group must be correct for the molecular geometry and normal modes. It is recommended to compare a low-quanta SA-VCI calculation with a full VCI calculation before using larger VCI spaces.

## ⚡ Performance Considerations

The size of the VCI space increases rapidly with the number of vibrational modes and the maximum number of quanta.

For large systems, consider using:

* A smaller `NQUANT`
* A limited number of states through `NSTATE`
* Selected VCI using `MAXSCI`
* Symmetry-adapted VCI using `PGROUP`
* Multiple OpenMP threads using `THREAD`

Full VCI calculations may require substantial memory because the Hamiltonian matrix is stored as a dense symmetric matrix.

## 🖱️ Graphical User Interface

The Windows GUI is available in the separate folder:

    GUI_and_precompiled_Windows/

The GUI can be used to:

* Select an ORCA `.vpt2` file
* Generate `input_vscf.txt`
* Run the precompiled Fortran executable
* Monitor calculation output in real time
* Stop running calculations
* Save calculation logs

The GUI and precompiled executable are intended for Windows only.

## 📊 Spectrum Viewer

The Windows folder also includes a graphical spectrum viewer that can read:

    vscf.out
    intensities.txt
    normal_mode.txt

The viewer can be used to:

* Plot harmonic, VSCF, and VCI spectra
* Apply Gaussian broadening
* Change spectral linewidth
* Apply temperature-dependent intensity corrections
* Load experimental JCAMP-DX spectra
* Perform baseline correction
* Inspect vibrational assignments
* Animate normal modes in three dimensions
* Export spectra as text files or PNG images

## 🔍 Troubleshooting

### The program cannot find the ORCA `.vpt2` file

Check the `FILECT` entry in `input_vscf.txt`. Use the complete path if the file is not located in the working directory.

### The calculation stops because of insufficient memory

Reduce `NQUANT`, reduce `NSTATE`, use Selected VCI through `MAXSCI`, or use symmetry adaptation through `PGROUP`.

### Symmetry-adapted VCI gives unexpected results

Verify that the selected point group matches the molecular structure. Run a small full VCI calculation and compare the energies with the SA-VCI results.

### The Windows executable does not run on Linux or macOS

The precompiled executable is Windows-only. Compile the Fortran source code on the target operating system.

## 📚 Citation

If you use ViBra in scientific work, please cite the associated publication and documentation.

## 📄 License

Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)

Copyright (c) 2026 Raphael F. Ligório and co-authors, CBPF (Centro Brasileiro de Pesquisas Físicas)

You are free to use, share, and modify ViBra for non-commercial purposes, provided you give appropriate credit. Commercial use is not permitted without explicit permission.

Full license text: http://creativecommons.org/licenses/by-nc/4.0/

## AI Usage

AI tools (LLMs e.g. Claude/DeepSeek) were used during development for code cleanup, formatting consistency, adding comments, and translation regarding bilingual (PT/EN) documentation. Where AI suggested boilerplate or repetitive code (e.g., expanding structurally similar routines), every such change was reviewed, tested, and verified by the developers before being incorporated. Core algorithms and numerical methods were designed, implemented, and validated independently of AI assistance.


## 📧 Contact

Remember: this is an under development project. For questions, bug reports, or contributions, please contact the project maintainers.
