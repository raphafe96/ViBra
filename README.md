# ViBra

ViBra is a program for anharmonic vibrational spectroscopy based on Vibrational Self-Consistent Field (VSCF), Vibrational Configuration Interaction (VCI), Selected VCI (S-VCI), and Symmetry-Adapted VCI (SA-VCI) methods.

The computational engine is written in Fortran 90/95 and reads anharmonic vibrational data from an ORCA `.vpt2` output file. ViBra can calculate vibrational energies, infrared intensities, normal modes, VSCF modals, VCI wavefunctions, and spectra.

## Main Features

* Harmonic oscillator, VSCF, VCI, Selected VCI, and Symmetry-Adapted VCI calculations
* Support for cubic and quartic force fields obtained from ORCA VPT2 calculations
* Infrared intensities using first- and second-order dipole derivatives
* OpenMP parallelization
* LAPACK/BLAS diagonalization routines
* Symmetry treatment for the Abelian point groups:

  * `C1`
  * `Cs`
  * `Ci`
  * `C2`
  * `C2h`
  * `C2v`
  * `D2`
  * `D2h`
* Output files compatible with the ViBra graphical spectrum viewer

## Repository Structure

```text
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
│   └── other Fortran source files
│
├── examples/
│   ├── input_vscf.txt
│   └── example ORCA .vpt2 files
│
├── GUI_and_precompiled_Windows/
│   ├── ViBra.exe
│   ├── GUI files
│   └── spectrum viewer files
│
└── README.md
```

## Important Note About the GUI and Precompiled Executable

The graphical user interface and the precompiled Fortran executable are provided in a separate folder:

```text
GUI_and_precompiled_Windows/
```

These files are intended to work only on Windows.

The precompiled executable was built for Windows, and the GUI was developed for use on Windows systems. Users working on Linux or macOS should compile the Fortran source code themselves and run the program from the command line.

## Requirements for Compilation

To compile ViBra from source, the following software is required:

* A Fortran compiler with Fortran 90/95 support
* OpenMP support
* LAPACK
* BLAS

The code was developed and tested using Intel Fortran (`ifx`). Other compilers may work, but the compilation flags and linked libraries may need to be adjusted.

Recommended compiler:

```text
Intel Fortran Compiler (ifx)
```

Recommended numerical libraries:

```text
Intel Math Kernel Library (MKL)
```

## Compilation on Windows Using Intel Fortran

An example compilation command using Intel Fortran is:

```bash
ifx /O3 /Qopenmp /threads /Qmkl:parallel /Qm64 /heap-arrays /fpscomp:logicals *.f90 /exe:ViBra.exe
```

This command enables optimization, OpenMP parallelization, threaded MKL routines, 64-bit compilation, and heap allocation for temporary arrays.

Depending on your installation, you may need to compile the source files in a specific order. A typical order is:

```text
integrals.f90
get_combination.f90
jacobi.f90
symmetry.f90
one_mode_operation.f90
read_input.f90
read_orca.f90
vci.f90
main.f90
```

If additional modules are included in the source directory, compile them before the files that use them.

## Compilation on Linux

A possible compilation command using Intel Fortran on Linux is:

```bash
ifx -O3 -qopenmp -qmkl=parallel -heap-arrays *.f90 -o ViBra
```

If using `gfortran`, LAPACK, BLAS, and OpenMP must be linked manually. For example:

```bash
gfortran -O3 -fopenmp *.f90 -llapack -lblas -o ViBra
```

The exact command may vary depending on the installed libraries and operating system.

## Input File

ViBra reads the input file:

```text
input_vscf.txt
```

The input file uses keyword-based entries. A typical example is:

```text
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
```

## Input Keywords

| Keyword  |    Type | Description                                                                                                          |
| -------- | ------: | -------------------------------------------------------------------------------------------------------------------- |
| `NMODES` | Integer | Number of vibrational modes.                                                                                         |
| `NEXPAN` | Integer | Number of harmonic oscillator basis functions used for each mode.                                                    |
| `FILECT` |  String | Path to the ORCA `.vpt2` output file.                                                                                |
| `CTEMOD` |  String | Input format. Use `orca_vpt2`.                                                                                       |
| `NQUANT` | Integer | Maximum total number of vibrational quanta included in the VCI space. Values less than or equal to zero disable VCI. |
| `NSTATE` | Integer | Number of eigenstates to calculate. Values less than or equal to zero calculate all states.                          |
| `CVGSCF` | Integer | VSCF convergence exponent. For example, `6` corresponds to a threshold of `10^-6 cm^-1`.                             |
| `THREAD` | Integer | Number of OpenMP threads.                                                                                            |
| `PGROUP` |  String | Molecular point group. Supported groups are `C1`, `Cs`, `Ci`, `C2`, `C2h`, `C2v`, `D2`, and `D2h`.                   |
| `PROJCT` |    Real | Projection cutoff used during symmetry analysis.                                                                     |
| `MAXSCI` | Integer | Number of configurations retained in Selected VCI. Use `0` for full VCI.                                             |

## Running ViBra

After compiling the code, place the executable in the same directory as `input_vscf.txt`, or provide the correct path to the input file and ORCA output file.

On Windows:

```bash
ViBra.exe
```

On Linux or macOS:

```bash
./ViBra
```

The program reads `input_vscf.txt` and starts the calculation.

## ORCA Input Requirement

ViBra requires an ORCA `.vpt2` output file containing the anharmonic vibrational information, including:

* Harmonic frequencies
* Normal modes
* Cubic force constants
* Quartic force constants
* First-order dipole derivatives
* Second-order dipole derivatives, when available

The ORCA calculation must be configured to generate the required VPT2 output.

## Output Files

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

## Calculation Modes

### Harmonic Oscillator Calculation

Set:

```text
NQUANT 0
```

This disables the VCI calculation and performs the harmonic and VSCF-related steps.

### Full VCI Calculation

Set:

```text
MAXSCI 0
```

and choose a positive value for `NQUANT`.

Example:

```text
NQUANT 4
MAXSCI 0
```

### Selected VCI Calculation

Set `MAXSCI` to a positive value.

Example:

```text
NQUANT 6
MAXSCI 100
```

A larger `MAXSCI` value retains more configurations and generally improves agreement with full VCI, at the cost of additional computational time.

### Symmetry-Adapted VCI Calculation

Set the molecular point group using `PGROUP`.

Example:

```text
PGROUP D2h
```

Symmetry-adapted VCI reduces computational cost by block-diagonalizing the VCI Hamiltonian according to irreducible representations.

The selected point group must be correct for the molecular geometry and normal modes. It is recommended to compare a low-quanta SA-VCI calculation with a full VCI calculation before using larger VCI spaces.

## Performance Considerations

The size of the VCI space increases rapidly with the number of vibrational modes and the maximum number of quanta.

For large systems, consider using:

* A smaller `NQUANT`
* A limited number of states through `NSTATE`
* Selected VCI using `MAXSCI`
* Symmetry-adapted VCI using `PGROUP`
* Multiple OpenMP threads using `THREAD`

Full VCI calculations may require substantial memory because the Hamiltonian matrix is stored as a dense symmetric matrix.

## Graphical User Interface

The Windows GUI is available in the separate folder:

```text
GUI_and_precompiled_Windows/
```

The GUI can be used to:

* Select an ORCA `.vpt2` file
* Generate `input_vscf.txt`
* Run the precompiled Fortran executable
* Monitor calculation output in real time
* Stop running calculations
* Save calculation logs

The GUI and precompiled executable are intended for Windows only.

## Spectrum Viewer

The Windows folder also includes a graphical spectrum viewer that can read:

```text
vscf.out
intensities.txt
normal_mode.txt
```

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

## Troubleshooting

### The program cannot find the ORCA `.vpt2` file

Check the `FILECT` entry in `input_vscf.txt`. Use the complete path if the file is not located in the working directory.

### The calculation stops because of insufficient memory

Reduce `NQUANT`, reduce `NSTATE`, use Selected VCI through `MAXSCI`, or use symmetry adaptation through `PGROUP`.

### Symmetry-adapted VCI gives unexpected results

Verify that the selected point group matches the molecular structure. Run a small full VCI calculation and compare the energies with the SA-VCI results.

### The Windows executable does not run on Linux or macOS

The precompiled executable is Windows-only. Compile the Fortran source code on the target operating system.

## Citation

If you use ViBra in scientific work, please cite the associated publication and documentation.

## License

Please add the appropriate license information for this project.

## Contact

For questions, bug reports, or contributions, please contact the project maintainers.
