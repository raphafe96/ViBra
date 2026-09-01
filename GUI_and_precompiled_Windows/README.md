## 🚀 Quick Start

1. Download `Final_dist.zip` from the repository.
2. Extract it to a desired location on your computer (do not use folder with spaces or special characters).
3. Download the precompiled Fortran engine `vscf_vci.exe`.
4. Place `vscf_vci.exe` **inside the extracted `Final_dist` folder**, in the same directory as `vscf_vci_gui.exe` and the other files.
5. Run `vscf_vci_gui.exe` to launch the main application.
6. It may take a while to launch the Python GUI interfaces, few seconds. Be patient here.

## 🖥️ Main Application (vscf_vci_gui.exe)

The main GUI provides input creation and execution of the Fortran engine. The interface is shown below:

![Main GUI - Vibra0](main/assets/Vibra0.png)

The GUI automatically detects the presence of `vscf_vci.exe` in the same folder and allows you to run calculations directly from the interface.

## 📊 Visualization (visualization.exe)

This separate executable is used for spectral and normal mode analysis. It offers two interfaces:

1. **Spectral viewer**  
   ![Spectral Viewer - Vibra1](https://raw.githubusercontent.com/raphafe96/ViBra/main/assets/Vibra1.png)

2. **Normal mode animator**  
   ![Normal Mode Animator - Vibra2](https://raw.githubusercontent.com/raphafe96/ViBra/main/assets/Vibra2.png)

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

### 🌐 Language Setting

If the GUI starts in Portuguese, it means the last language used was Portuguese. To switch to English, in the vscf_vci_gui.exe click the purple button labeled **"Idioma"** (which means "Language" in Portuguese). This toggles the interface language.

