#!/usr/bin/env bash
# Set up the PyVista plotting environment for PyClical plotting demos.
#
# Run from the glucat repository root with:
#   source pyclical/demos/plotting/setup-pyvista-env.sh
#
# On x86-64:
#   Creates/updates the 'pyclical-pyvista' Conda environment from
#   pyclical/demos/plotting/pyvista-env.yml and activates it.
#
# On ARM aarch64 (Asahi Linux):
#   Creates a --system-site-packages virtual environment at .venvs/pyclical-pyvista,
#   activates it, installs PyVista without bundled VTK (pip install --no-deps pyvista)
#   to utilize the 16KB page-aligned system python3-vtk RPM, installs auxiliary
#   requirements from requirements-pyvista.txt, and configures the lib64 search path.

if [ "${BASH_SOURCE-}" = "$0" ]; then
    echo "Please run this script from Bash as follows:" >&2
    echo "  source $0" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.."; pwd)"
_ARCH=$(uname -m)

if [ "$_ARCH" = "aarch64" ]; then
    echo "Configuring PyVista environment for ARM64 (aarch64) system site-packages ..."
    if ! python3 -c 'import vtkmodules' 2>/dev/null; then
        echo "Error: python3-vtk is not installed in system site-packages." >&2
        echo "The 16KB page-aligned system VTK RPM is required on ARM64 / Asahi Linux." >&2
        echo "Please run: sudo dnf install python3-vtk" >&2
        return 1 2>/dev/null || exit 1
    fi

    mkdir -p "${REPO_ROOT}/.venvs"
    python3 -m venv --system-site-packages "${REPO_ROOT}/.venvs/pyclical-pyvista"
    source "${REPO_ROOT}/.venvs/pyclical-pyvista/bin/activate"

    echo "Installing PyVista without bundled VTK (reusing system python3-vtk RPM)..."
    pip install --no-deps pyvista
    pip install -r "${REPO_ROOT}/pyclical/demos/plotting/requirements-pyvista.txt"

    # Inject lib64 path into venv purelib if needed (Fedora lib64 site-packages path)
    python3 -c "
import sysconfig, os, sys
purelib = sysconfig.get_path('purelib')
if purelib and os.path.exists(purelib):
    pth = os.path.join(purelib, 'fedora_lib64.pth')
    lib64 = '/usr/lib64/python{}.{}/site-packages'.format(*sys.version_info[:2])
    with open(pth, 'w') as f:
        f.write(lib64 + '\n')
    print('Configured venv site-packages path:', pth, '->', lib64)
"
    echo ""
    echo "Virtual environment '.venvs/pyclical-pyvista' is ready."
else
    if ! command -v mamba >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
        echo "Error: neither mamba nor conda found on PATH." >&2
        echo "Install Miniforge (https://github.com/conda-forge/miniforge) and try again." >&2
        return 1
    fi

    if command -v mamba >/dev/null 2>&1; then
        _CONDACMD=mamba
    else
        _CONDACMD=conda
    fi

    echo "Creating/updating Conda environment from pyclical/demos/plotting/pyvista-env.yml ..."
    $_CONDACMD env create -f "${REPO_ROOT}/pyclical/demos/plotting/pyvista-env.yml" 2>/dev/null \
        || $_CONDACMD env update -f "${REPO_ROOT}/pyclical/demos/plotting/pyvista-env.yml"

    $_CONDACMD activate pyclical-pyvista

    echo ""
    echo "Environment 'pyclical-pyvista' is ready."
fi

echo "Next steps (from the repository root, with this environment active):"
echo "  make -f admin/Makefile.common bootstrap  # git clone only, not needed for tarballs"
echo "  ./configure"
echo "  make -C pyclical -j$(( $(nproc)/2 ))"
echo "  source pyclical/demos/plotting/export-pyvista-vars.sh"
