# 🔬 Welcome to the ViBra Playground!

This is the ViBra playground, a sandbox for new, experimental, and evolving features that push the boundaries of vibrational configuration interaction. Explore, experiment, and help shape our software. Just remember: things here may change, break, or evolve, that is the nature of a playground!

## ⚠️ Important Warning: Selected VCI with List Mode (State List Parsing Bug)

A severe bug was found when the Selected VCI is used with **list mode** (`MAXSCI N list`). The issue arises in the construction of the CI reference state: ViBra builds the list of total configurations based on total quanta, then reads the user‑provided state list and compares it with all states to define what is reference and what is external (for further EN-PT enlargement, if requested). The bug was in the ordering/labeling of the states: the program was using the order of the full configuration list rather than the order of the provided list. As a result, the reference space contained **the exact same number** of states as in the list file, but because the ordering was not properly parsed, the final reference contained **mismatched states** from the actual list provided.

This does not affect the overall trend of the results presented in the arXiv paper and Zenodo repository, but it does affect the exact numerical values. We are waiting for further validation before updating the arXiv paper and Zenodo repository. A working (but still being tested) version is available here.

✅ **Note:** The **auto mode** for Selected VCI (`MAXSCI N auto s/d/t/q`) is **working as intended** and is **not affected** by the list‑mode parsing bug.

Update: the Zenodo repository and the source code have been replaced with a working version/correct results.

## ✨ What's New in This Version

The playground currently introduces the following new and experimental features:

- **VCI@HO**: Possibility to run VCI using a harmonic oscillator basis for all VCI routines, including Full VCI, Selected VCI (S-VCI), and Symmetry-Adapted VCI (SA-VCI). This provides an alternative to the standard VCI@VSCF approach and can serve as a useful reference.
- **Iterative diagonalization**: Option to run an iterative Davidson diagonalizer when a full VCI or selected VCI calculation is requested, avoiding the memory bottleneck of storing the full dense Hamiltonian matrix. Not yet implemented for the symmetry adapted part.
- **Built-in water molecule test**: An option to automatically run a water molecule test case using VCI@VSCF. This allows you to compare results with Crystal23 output and verify that everything is running correctly (you must explicitly set the RUNSCF to 1).
- **Extended intensity output**: For full and selected VCI calculations, frequencies are now also saved in km/mol alongside the transition dipoles in a new file called `dipoles_intensity_vci.txt`, providing richer data for spectral analysis.
- **Mode exclusion**: Possibility to exclude specific vibrational modes from the VSCF/VCI/VPT calculation entirely, either automatically (below a frequency cutoff) or by explicitly listing mode indices. **This currently does not work for Symmetry-Adapted VCI (SA-VCI)** — see the keyword description and warning below.
- **Sum over States VPT2 (SoS VPT2)**: A way to obtain anharmonic energies without a full VCI diagonalization, using the existing VCI Hamiltonian kernel. Described in detail below.
- **Force-field term exclusion by index-distinctness**: Possibility to zero out specific cubic/quartic force constants before VSCF/VCI/VPT2, based on how many distinct mode indices they involve. Described in detail below.

## 🧪 How to Use the New Features

To activate the playground features, simply place a file named `extra_input.txt` in the same folder as your `input_vscf.txt`. The keywords in this file will be merged into the normal input. Once these features graduate to the stable release, they will be integrated directly into `input_vscf.txt`.

### Extra Input Keywords

| Keyword  |    Type | Description                                                                                                         |
| -------- | ------: | ------------------------------------------------------------------------------------------------------------------- |
| `RUNSCF` | Integer | Set to 1 to run as VCI@VSCF, or 0 to run as VCI@HO. Default: 1.                                                               |
| `RUNH2O` | Integer | Set to 1 to run a test mode for the water molecule, or 0 for normal runs. Default: 0.                                          |
| `RUNPT2` | Integer | Set to 1 to run Sum over States VPT2 (SoS VPT2) instead of a VCI diagonalization, or 0 for a normal VCI run. Default: 0. See the dedicated section below for what this changes and why it exists. |
| `RUNDAV` | Integer | Set to 1 to use the iterative Davidson diagonalizer, or 0 to use standard dense matrix diagonalization. Default: 0         |
| `DAVCUT` |    Real | Cutoff threshold for Hamiltonian matrix elements in the Davidson solver. Default: 0.0000005.                       |
| `DAVCVG` |    Real | Convergence threshold for the iterative diagonalizer. Default: 0.0001.                                              |
| `DAVMAX` | Integer | Maximum number of iterations for the Davidson solver. Default: 25.                                                 |
| `DAVSTA` | Integer | Number of lowest eigenvalues to compute. Default: 10.                                                               |
| `RUNENR` | Integer | Set to 1 to estimate the number of states needed for a given frequency threshold (based on HO energies), 2 to also print the state energies (HO), or 0 to do nothing. Default: 0. |
| `MAXFRQ` |    Real | Maximum frequency in cm⁻¹ for which intensities will be calculated. Default: 4500.0.                                |
| `DAVBUF` | Integer | Buffer size for the subspace dimension in the Davidson diagonalizer. There is a minimum limit internally set to DAVSTA times 20. Default: 4000. |
| `EXCLUD` |  Mixed  | Excludes vibrational modes from the calculation before VSCF/VCI. Two sub-keyword forms: `EXCLUD auto <freq_cutoff>` removes every mode with a harmonic frequency (cm⁻¹) below `<freq_cutoff>` (real); `EXCLUD spec <mode1> <mode2> ...` (integers) removes exactly the listed mode indices. For `spec`, indices refer to **vibrational** modes only: a non-linear molecule with N atoms has 3N total modes, of which 3N−6 are vibrational after removing the 3 translations and 3 rotations, and the first vibrational mode is index 1 (not the first of the 3N raw modes). Default: not set (no exclusion). **⚠️ Does not currently work with Symmetry-Adapted VCI (SA-VCI, i.e. a point group other than C1).** |
| `R3DIFF` | Integer | Set to 1 to zero every cubic force constant Φ_ijk with 3 distinct mode indices (i≠j≠k≠i), before VSCF/VCI/VPT2. Default: 0. |
| `R4DIFF` | Integer | Set to 1 to zero every quartic force constant Φ_ijkl with 4 distinct mode indices, before VSCF/VCI/VPT2. Default: 0. |
| `R4TRIP` | Integer | Set to 1 to zero every quartic force constant Φ_ijkl with exactly 3 distinct mode indices (the Φ_iijk-type patterns), before VSCF/VCI/VPT2. Default: 0. |

### Sum over States VPT2 (SoS VPT2)
See: https://dx.doi.org/10.1021/acs.jpca.0c09526, J. Phys. Chem. A 2021, 125, 1301−1324'

Setting `RUNPT2 1` switches the dispatch away from VCI entirely and calls a dedicated `vpt2` routine that computes second order perturbative energies directly, using a sum over states formulation:

    E(i) = U(i) + <Ψ(i)|V|Ψ(i)> + Σ_{j≠i} |<Ψ(i)|V|Ψ(j)>|² / (U(i) − U(j))

where `U(i)` is the zero order harmonic energy of state `i` and `V` is the anharmonic (cubic and quartic) part of the potential.

The motivation for this route is practical rather than purely theoretical: the matrix element kernel used to build the sparse VCI Hamiltonian (the routine that evaluates `<Ψm|Ĥ|Ψn>` from the cubic and quartic force constants and the precomputed modal integrals) was already implemented and tested for VCI. SoS VPT2 reuses that same kernel unchanged: it builds the identical sparse pair list of nonzero couplings, evaluates the same matrix elements, and then, instead of diagonalizing them, plugs the diagonal elements and the off diagonal couplings directly into the second order perturbative expression above. No new integral machinery was needed to add this feature.

Because this formulation assumes the basis states diagonalize the zero order (harmonic) Hamiltonian exactly, the code automatically switches the ground state modals to the harmonic oscillator basis whenever `RUNPT2` is nonzero, overriding whatever value `RUNSCF` was given in the input. In other words, requesting SoS VPT2 always runs on top of VCI@HO style modals, never VCI@VSCF ones, and you do not need to set `RUNSCF 0` yourself, it happens automatically.

This sum over states construction differs from the standard, formally ordered VPT2 expression used by packages (which keeps only the diagonal quartic contribution at first order and the off diagonal cubic contribution at second order, dropping off diagonal quartic terms as higher order in the perturbation parameter). SoS VPT2 as implemented here keeps the off diagonal quartic contributions as well, so it is best thought of as a related but distinct second order treatment rather than a drop in replacement for standard VPT2 output.

### Force-Field Term Exclusion by Index-Distinctness

Force fields used here are truncated at fourth order and contain all cubic and quartic terms up to 4 distinct mode indices. Some other codes and workflows instead build (or are limited to) a reduced representation of the PES that neglects some or all of the higher-order mode-coupling terms. Neglecting these couplings almost always results in a wrong description of the vibrational states, so these keywords are provided primarily as a **diagnostic and comparison tool**, not as a recommended production setting.

Three independent keywords control which coupling terms are zeroed, based purely on how many distinct mode indices a given force constant involves:

* `R3DIFF 1` — zeros cubic terms Φ_ijk with all 3 indices distinct.
* `R4DIFF 1` — zeros quartic terms Φ_ijkl with all 4 indices distinct.
* `R4TRIP 1` — zeros quartic terms Φ_ijkl with exactly 3 distinct indices (the Φ_iijk-type patterns).

They can be set independently or in any combination, e.g. setting all three together drops every cubic/quartic term with 3 or more distinct mode indices, keeping only terms with at most 2 distinct indices (Φ_iii, Φ_iiii, Φ_iij, Φ_iiij, Φ_ijjj, Φ_iijj).

The exclusion is applied to `Potential_3`/`Potential_4` right after the degeneracy factors are applied and before the sparse inverted index is built, so any dropped terms are fully absent from every downstream routine — VSCF, VCI, and SoS VPT2 alike (`RUNPT2 1` reuses the same reduced `Potential_3`/`Potential_4` arrays). The number of zeroed cubic and quartic entries is reported to both stdout and `vscf.out`.

**Note on `R4DIFF`:** ORCA performs a semi-quartic force field, where fully 4-distinct-index terms (Φ_ijkl) don't exist, and the ORCA reader doesn't even parse anything of that form. `R4DIFF` is kept for the future, in case force fields from other software with a genuine full quartic representation are supported later.

**Why remove real coupling terms at all?** 

* **Diagnosing SCF instability.** If VSCF diverges or oscillates, selectively zeroing term classes helps isolate which coupling is responsible (e.g. an unbounded or ill-conditioned Φ_iijj), before deciding whether the fix belongs in the force field itself or elsewhere.
* **Numerical unreliability of mixed derivatives.** Force constants mixing 3+ distinct modes come from mixed finite-difference derivatives, which accumulate more numerical noise than pure single-mode (diagonal) derivatives — smaller signal, more cancellation error, more sensitivity to displacement step size. If a specific mixed term is noise-dominated rather than physically meaningful, removing it discards numerical error rather than real physics.
* **Matching another method or code for benchmarking.** Some codes or reference datasets only ever compute a reduced-coupling force field; truncating your own the same way gives a fair, like-for-like comparison rather than crediting your calculation with physics the reference never had.
* **Cost at large mode counts.** The number of 3-/4-distinct terms scales combinatorially with `N_modes`; for large systems, dropping a subset that is small and/or unreliable can be a pragmatic trade-off alongside the physical justification above.

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

### Example: Running SoS VPT2

To run Sum over States VPT2 instead of a VCI diagonalization, create an `extra_input.txt` file with:

    RUNPT2 1

Configure `input_vscf.txt` as usual (NQUANT, force field files, and so on); the ground state modals will automatically be the harmonic oscillator ones regardless of the `RUNSCF` value.

### Example: Excluding Modes

To automatically exclude all modes below 200 cm⁻¹:

    EXCLUD auto 200.0

To exclude specific modes (e.g. modes 1, 6, and 38):

    EXCLUD spec 1 6 38

⚠️ Only use `EXCLUD` with point group `C1` for now. Combining it with SA-VCI (any other point group) will give incorrect irrep assignments and, consequently, incorrect symmetry-restricted VCI results.

### Example: 0 Force-Field Terms by Index-Distinctness

To zero only the fully off-diagonal cubic terms (Φ_ijk, 3 distinct indices):

    R3DIFF 1

To drop all cubic and quartic terms with 3 or more distinct mode indices, keeping only terms with at most 2 distinct indices:

    R3DIFF 1
    R4DIFF 1
    R4TRIP 1

This is primarily useful as a diagnostic: comparing full-force-field results against a reduced-term run isolates how much of a given spectral feature comes from genuine higher-mode coupling versus the lower-order terms. It is not a substitute for a converged calculation with the full force field.

## Changelog

**29/08/2026**

* Added `found_hessian`, `found_cubic`, `found_quartic`, `found_dipole1`, and `found_dipole2` logical flags to the `read_orca` subroutine (`read_orca.f90`), set as each corresponding block is located while parsing `basename.vpt2`. If any expected block (Hessian, cubic force constants, quartic force constants, first dipole derivatives, or second dipole derivatives) is missing from the file, the program now stops with a specific error message naming the exact header line that was not found, instead of failing later with an unrelated or unclear error.
* Added the `R3DIFF`, `R4DIFF`, and `R4TRIP` keywords to `extra_input.txt` and a new `zero_offdiagonal_mode_terms` subroutine in `get_combination.f90`: each independently zeros a specific class of cubic/quartic force constants (Φ_ijk with 3 distinct indices, Φ_ijkl with 4 distinct indices, and Φ_ijkl with exactly 3 distinct indices, respectively) before VSCF/VCI/VPT2. Reports the number of zeroed cubic and quartic entries to stdout and `vscf.out`.

**27/08/2026**

* Fixed a bug in the parsing of the user‑provided state list in the `selected_vibrational_ci` subroutine (`vci.f90`), specifically in step 3 (building the CI reference state). The previous version incorrectly used the ordering of the full configuration list instead of the user‑provided list, causing mismatched reference states.
* Added a checker to detect and reject duplicated states in the user‑provided list.
* The keyword `NSTATE` is now **deprecated**. If fewer than the total number of states is desired, run a iterative diagonalization and set the number of roots (`DAVSTA`) accordingly.

**25/08/2026**

* Added a new file `vpt2.f90`, containing the new `v_pt2` module and its `vpt2` subroutine, and added the new `RUNPT2` keyword in `main_vscf.f90` to dispatch to it, implementing Sum over States VPT2 (SoS VPT2):

  * Reuses the existing VCI sparse Hamiltonian kernel to evaluate diagonal and off diagonal matrix elements, then applies second order perturbation theory instead of diagonalizing.
  * Automatically forces the harmonic oscillator ground state modals (equivalent to `RUNSCF 0`) whenever `RUNPT2` is nonzero.

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
