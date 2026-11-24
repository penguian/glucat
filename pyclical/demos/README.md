# PyClical Demos and Tutorials

This directory contains demonstration scripts, interactive tutorials, Jupyter notebooks,
Marimo reactive notebooks, and 3D visualization applications demonstrating the capabilities
of PyClical (the Python interface to GluCat).

## Directory Organization

The demos are divided into two focused subdirectories:

- **`sandbox/`**: Standard tutorials, presentation scripts, reference outputs, and computational notebooks.
  - `pyclical_tutorials.py`: Interactive CLI runner providing a menu for tutorials (`0_0` through `2_0`).
  - `clifford_demo.py`: Clifford algebra presentation demo with offline portrait rendering.
  - `m_theory_demo.py`: 11-dimensional M-theory Clifford algebra demonstration.
  - `pyclical_demo.py` & `sqrt_log_demo.py`: Core functionality and matrix logarithm/square root benchmarks.
  - `versor_demo.py`: Versor and rotor demonstrations.
  - `build_pyclical_notebooks.py`: Builds Jupyter notebooks (`.ipynb`) from tutorial scripts.
  - `build_marimo_notebooks.py`: Builds reactive Marimo notebooks (`*_marimo.py`) from tutorial scripts.
  - `pyclical_marimo_utils.py`: Runtime helpers for Marimo notebooks (sequential PyClical session execution, clickable `▶` output toggles, and interactive exercise grading).
  - `validate_notebooks.py`: Automated validation test runner for generated Jupyter notebooks.
  - `images/`: Offline historical portraits used by notebooks and presentation scripts.
- **`plotting/`**: 3D geometric visualization using PyVista and PySide6.
  - `plotting_demo_pyvista.py`: 3D orbits and torus visualization using PyVista.
  - `plotting_demo_pyvista_dialog.py`: Interactive PySide6 GUI window with parameter sliders for PyVista visualization.
  - `plotting_demo_dialog.py`: Historical reference dialog (BSD licensed).

## Quick Start

### 1. Interactive Tutorials (CLI)
From the repository root or sandbox directory:
```bash
python3 pyclical/demos/sandbox/pyclical_tutorials.py
```
This presents an interactive menu to select and execute any tutorial from 0.0 to 2.0.

### 2. Presentation Scripts
```bash
cd pyclical/demos/sandbox
python3 clifford_demo.py
python3 m_theory_demo.py
python3 versor_demo.py
```

### 3. JupyterLab Notebooks
To explore and execute the tutorials inside JupyterLab:
```bash
jupyter lab pyclical/demos/sandbox/
```
*(If using system Python, register the PyClical kernel first via `make -C pyclical install-pyclical-kernel`).*

### 4. Marimo Reactive Notebooks
Ensure `marimo` is installed in the active Python environment used to compile PyClical (included automatically by `source pyclical/setup-env.sh`, or install via `conda install -c conda-forge marimo` or `pip install marimo`). Unlike Jupyter, Marimo does not use `kernel.json` kernel specs; it executes directly inside the Python environment that runs the `marimo` CLI.

Each generated notebook displays immutable syntax-highlighted code blocks with clickable `▶` output toggles, plus interactive `Check` exercise editors and collapsible `"Show solution"` accordions:
```bash
# Generate reactive Marimo notebooks (*_marimo.py):
python3 pyclical/demos/sandbox/build_marimo_notebooks.py

# Run in read-only interactive App mode:
marimo run pyclical/demos/sandbox/clifford_demo_marimo.py

# Or edit and inspect cells interactively in Studio mode:
marimo edit pyclical/demos/sandbox/clifford_demo_marimo.py
```

### 5. 3D Plotting Demos (PyVista)
PyVista requires graphics rendering libraries and Qt bindings:
```bash
# Activate environment (or ensure pyvista and pyside6 are installed)
source pyclical/setup-env.sh

cd pyclical/demos/plotting
python3 plotting_demo_pyvista.py
python3 plotting_demo_pyvista_dialog.py

# To run headless/non-interactive test:
GLUCAT_NON_INTERACTIVE=1 python3 plotting_demo_pyvista.py
```
