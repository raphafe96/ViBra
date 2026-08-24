# 🔬 Welcome to the ViBra Playground!

This is the ViBra playground, a sandbox for new, experimental, and evolving features that push the boundaries of vibrational configuration interaction. Explore, experiment, and help shape our software. Just remember: things here may change, break, or evolve, that is the nature of a playground!

## ✨ What's New in This Version

The playground currently introduces the following new and experimental features:

- **VCI@HO**: Possibility to run VCI using a harmonic oscillator basis for all VCI routines, including Full VCI, Selected VCI (S-VCI), and Symmetry-Adapted VCI (SA-VCI). This provides an alternative to the standard VCI@VSCF approach and can serve as a useful reference.
- **Iterative diagonalization**: Option to run an iterative Davidson diagonalizer when a full VCI or selected VCI calculation is requested, avoiding the memory bottleneck of storing the full dense Hamiltonian matrix. Not yet implemented for the symmetry adapted part.
- **Built-in water molecule test**: An option to automatically run a water molecule test case using VCI@VSCF. This allows you to compare results with Crystal23 output and verify that everything is running correctly (you must explicitly set the RUNSCF to 1).
- **Extended intensity output**: For full and selected VCI calculations, frequencies are now also saved in km/mol alongside the transition dipoles in a new file called `dipoles_intensity_vci.txt`, providing richer data for spectral analysis.
- **Mode exclusion (`EXCLUD`)**: Possibility to exclude specific vibrational modes from the VSCF/VCI calculation entirely, either automatically (below a frequency cutoff) or by explicitly listing mode indices. **This currently does not work for Symmetry-Adapted VCI (SA-VCI)** — see the keyword description and warning below.

## 🧪 How to Use the New Features

To activate the playground features, simply place a file named `extra_input.txt` in the same folder as your `input_vscf.txt`. The keywords in this file will be merged into the normal input. Once these features graduate to the stable release, they will be integrated directly into `input_vscf.txt`.

### Extra Input Keywords

| Keyword  |    Type | Description                                                                                                         |
| -------- | ------: | ------------------------------------------------------------------------------------------------------------------- |
| `RUNSCF` | Integer | Set to 1 to run as VCI@VSCF, or 0 to run as VCI@HO. Default: 1.                                                               |
| `RUNH2O` | Integer | Set to 1 to run a test mode for the water molecule, or 0 for normal runs. Default: 0.                                          |
| `RUNDAV` | Integer | Set to 1 to use the iterative Davidson diagonalizer, or 0 to use standard dense matrix diagonalization. Default: 0         |
| `DAVCUT` |    Real | Cutoff threshold for Hamiltonian matrix elements in the Davidson solver. Default: 0.0000005.                       |
| `DAVCVG` |    Real | Convergence threshold for the iterative diagonalizer. Default: 0.0001.                                              |
| `DAVMAX` | Integer | Maximum number of iterations for the Davidson solver. Default: 25.                                                 |
| `DAVSTA` | Integer | Number of lowest eigenvalues to compute. Default: 10.                                                               |
| `RUNENR` | Integer | Set to 1 to estimate the number of states needed for a given frequency threshold (based on HO energies), 2 to also print the state energies (HO), or 0 to do nothing. Default: 0. |
| `MAXFRQ` |    Real | Maximum frequency in cm⁻¹ for which intensities will be calculated. Default: 4500.0.                                |
| `DAVBUF` | Integer | Buffer size for the subspace dimension in the Davidson diagonalizer. There is a minimum limit internally set to DAVSTA times 20. Default: 4000. |
| `EXCLUD` |  Mixed  | Excludes vibrational modes from the calculation before VSCF/VCI. Two sub-keyword forms: `EXCLUD auto <freq_cutoff>` removes every mode with a harmonic frequency (cm⁻¹) below `<freq_cutoff>` (real); `EXCLUD spec <mode1> <mode2> ...` (integers) removes exactly the listed mode indices. For `spec`, indices refer to **vibrational** modes only: a non-linear molecule with N atoms has 3N total modes, of which 3N−6 are vibrational after removing the 3 translations and 3 rotations, and the first vibrational mode is index 1 (not the first of the 3N raw modes). Default: not set (no exclusion). **⚠️ Does not currently work with Symmetry-Adapted VCI (SA-VCI, i.e. a point group other than C1).** |

### Example: Running a Water Test

To run the built-in water molecule test and verify your setup, create an `extra_input.txt` file with the following content:

    RUNSCF 1
    RUNH2O 1

Then run ViBra as usual. The program will automatically execute a VCI@VSCF calculation for water. Compare the output energies and intensities with your Crystal23 reference to confirm everything is functioning properly.

### Example: Exploring VCI@HO with On-the-fly Diagonalization

To explore a large VCI space using a harmonic basis and the memory-efficient Davidson solver, create an `extra_input.txt` file with:

    RUNSCF 0
    RUNDAV 1
    DAVSTA 50

Then configure your standard `input_vscf.txt` with a large NQUANT and MAXSCI set to 0 for a full VCI. Run ViBra and the calculation will proceed in the VCI@HO framework, diagonalizing the Hamiltonian iteratively and saving significant memory.

### Example: Excluding Modes

To automatically exclude all modes below 200 cm⁻¹:

    EXCLUD auto 200.0

To exclude specific modes (e.g. modes 1, 6, and 38):

    EXCLUD spec 1 6 38

⚠️ Only use `EXCLUD` with point group `C1` for now. Combining it with SA-VCI (any other point group) will give incorrect irrep assignments and, consequently, incorrect symmetry-restricted VCI results.


## Changelog

**24/08/2026**

* Added the new `read_exclude` subroutine to `read_input.f90` and updated `main_vscf.f90`:

  * Added the `EXCLUD` keyword to `extra_input.txt`, allowing modes to be excluded from the VSCF/VCI calculation automatically (`auto`, by frequency cutoff) or explicitly (`spec`, by vibrational mode index).
  * Not yet compatible with Symmetry-Adapted VCI (`symmetry.f90` does not account for excluded modes).

**19/08/2026**

* Updated the `selected_vibrational_ci` subroutine in `vci.f90`:

  * Removed obsolete and unused variables.
  * Added intensities in km/mol.
  * Implemented iterative diagonalization.
  * Improved sparse pair-list construction for ORCA semi-quartic force fields (`ndif <= 3`).

**18/08/2026**

* Replaced `qsort_diag` with `qsort_diag2` in `jacobi.f90` (introduced a new sorting algorithm).

**18/08/2026**

* Started documenting the changelog here. For previous changes, check the code backups.

**July–August 2026**

* Added iterative diagonalization for full VCI.
* Added intensities in km/mol for full VCI.
* Added a simple water-molecule test case.
* Added VCI@HO.
