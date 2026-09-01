## 🚀 Quick Start

1. Download `Final_dist.zip` from the repository.
2. Extract it to a desired location on your computer (do not use folder with spaces or special characters).
3. Download the precompiled Fortran engine `vscf_vci.exe`.
4. Place `vscf_vci.exe` **inside the extracted `Final_dist` folder**, in the same directory as `vscf_vci_gui.exe` and the other files.
5. Run `vscf_vci_gui.exe` to launch the main application.
6. It may take a while to launch the Python GUI interfaces, few seconds. Be patient here.

## 🖥️ Main Application (vscf_vci_gui.exe)

The main GUI provides input creation and execution of the Fortran engine. The interface is shown below:

<img src="https://raw.githubusercontent.com/raphafe96/ViBra/main/assets/Vibra0.png" alt="Main GUI - Vibra0" width="500"/>

The GUI automatically detects the presence of `vscf_vci.exe` in the same folder and allows you to run calculations directly from the interface.

### 🌐 Language Setting

If the GUI starts in Portuguese, it means the last language used was Portuguese. To switch to English, in the `vscf_vci_gui.exe` click the purple button labeled **"Idioma"** (which means "Language" in Portuguese) at the top right corner. This toggles the interface language.

### 🧩 Graphical User Interface (Keyword Input)

The GUI simplifies job preparation. The user navigates to the ORCA `.vpt2` file; the working directory is set automatically. A modal dialog collects all keywords, validates their ranges with tooltip guidance, and writes `input_vscf.txt`. The Fortran executable is launched as a subprocess in a background thread, with real-time output streaming to a scrollable panel. A stop button terminates a running job, and the panel supports save and clear operations.

#### Principal Input Keywords

| Keyword | Type | Description |
|---------|------|-------------|
| `NMODES` | int | Number of vibrational modes (M). |
| `NEXPAN` | int | HO basis size per mode (N_exp). |
| `FILECT` | str | Path to ORCA `.vpt2` file. |
| `CTEMOD` | str | Format (`orca_vpt2`). |
| `NQUANT` | int | Maximum total quanta (N_q); ≤ 0 disables VCI. |
| `NSTATE` | int | Number of eigenstates (≤ 0 = all). |
| `CVGSCF` | int | VSCF convergence exponent (10^(-C) cm⁻¹). |
| `THREAD` | int | Number of OpenMP threads. |
| `PGROUP` | str | Point group (C1, Cs, Ci, C2, C2h, C2v, D2, D2h). |
| `PROJCT` | real | Projection cutoff for symmetry detection (Å). |
| `MAXSCI` | int | N_sel for S-VCI. Default: 100.<br>Usage: `MAXSCI [N] [mode] [ref]`<br>• N > 0, mode = auto: iterative selection with CI reference (ref = s, d, t, q).<br>• N > 0, mode = list: iterative selection with list reference, no ref needed.<br>• N = 0: full VCI (no S-VCI), no mode and no ref needed.<br>• N = 0, mode = list: full VCI using user-provided state list, no ref needed. |

### 📁 Output Files

The Fortran engine produces three output files, which can be utilized by the visualisation interface:

- `vscf.out` — complete log: VSCF convergence histories, modal coefficient matrices, VCI configuration list, eigenstate energies, and the three leading CI coefficients with quantum-number assignments;
- `intensities.txt` — transition frequencies (cm⁻¹) and normalised intensities for HO, VSCF, and VCI/SCI methods;
- `normal_mode.txt` — equilibrium geometry (Å) and Cartesian displacement vectors for all 3N modes.

## 📊 Visualization (visualization.exe)

This separate executable is used for spectral viewer:

<img src="https://raw.githubusercontent.com/raphafe96/ViBra/main/assets/Vibra1.png" alt="Main GUI - Vibra1" width="500"/>

### 📈 Interactive Spectral Viewer

The viewer (`MainApplication`, Python 3, Matplotlib, 3Dmol.js) reads the three output files and provides a rich environment for analysis and comparison:

- **Gaussian broadening:** each transition at frequency νₖ with intensity Iₖ is represented as  
  S(ν) = Σₖ Iₖ exp[−(ν−νₖ)²/(2σ²)],  
  with σ = FWHM/(2√(2 ln 2)). FWHM is adjustable via stepper and direct text entry.
- **Temperature-dependent spectra:** for T > 0 K, VCI intensities are Boltzmann-weighted:  
  Iₖᵉᶠᶠ = Iₖ (p₀ − pₖ), where pₖ = exp(−h c νₖ / k_B T)/Z and hc/k_B = 1.438777 cm·K.  
  Temperature is adjustable from 0 to 950 K.
- **JCAMP-DX overlay:** experimental spectra in JCAMP-DX format (`.jdx`/`.dx`) can be loaded. Transmittance is automatically converted to absorbance, negative values are clipped, and the result is normalised.
- **ALS baseline correction:** for experimental spectra exhibiting a sloping baseline, an Asymmetric Least Squares (ALS) algorithm estimates and subtracts a smooth baseline. The cost function  
  Q(z) = Σᵢ wᵢ(yᵢ − zᵢ)² + λ‖Dz‖²  
  is minimised iteratively, where D is the second-difference matrix. The smoothness parameter λ and asymmetry parameter p are adjustable via steppers with real-time redisplay.
- **Peak inspection:** clicking a VCI peak label displays the frequency, intensity, and three leading CI coefficients with quantum-number assignments (e.g., "v₁ + 2v₃"). Buttons launch 3D normal-mode animations for any involved mode.
- **Data export:** the displayed spectrum can be saved as a text file and as a 300 dpi PNG image.
- **3D normal-mode animation:** sinusoidal displacement trajectories (60 frames, adjustable amplitude) are generated from `normal_mode.txt`, saved as JSON, and rendered via an auto-generated HTML page with embedded 3Dmol.js viewer supporting play/pause, frame slider, speed control, and adjustable rendering styles and atomic radii. If a symmetry-adapted VCI is performed, the symmetry elements can also be displayed. 


  <img src="https://raw.githubusercontent.com/raphafe96/ViBra/main/assets/Vibra3.png" alt="Main GUI - Vibra3" width="500"/>

## ⚠️ Important Notes

### ⚙️ Precompiled Fortran Engine (vscf_vci.exe)

- The Fortran engine is distributed separately because it can be updated frequently.
- **Always check** that you have the latest version of `vscf_vci.exe` by consulting the repository.
- You may also compile the engine yourself, especially if you wish to use new features present in the [`playground`](https://github.com/raphafe96/ViBra/tree/main/playground) folder.
- Note: New input keywords introduced in the `playground` or source code may not yet be available in the GUI input generator. They will be added in future GUI versions.

### 🔧 Compilation of vscf_vci.exe

The precompiled engine is built from the source code in [`source`](https://github.com/raphafe96/ViBra/tree/main/source). It uses the LAPACK/BLAS wrapper from Intel MKL and is compiled with the Intel Fortran compiler `ifx`, using heap arrays and `-O3` optimization.

### 📖 Always check the main folder for the most recent manual

The GUI includes a button to open the manual (available in both Portuguese and English). However, the manual version included in the GUI may be outdated because it was bundled at the time the GUI was last compiled. **Always check the main repository folder for the most recent version of the manual.**
