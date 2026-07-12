#!/usr/bin/env bash
# Set up the pyclical-plotting Conda environment for Mayavi plotting demos.
#
# Run from the glucat repository root with:
#   source pyclical/demos/plotting/setup-plotting-env.sh
#
# This script:
#   1. Creates (or updates) the pyclical-plotting Conda environment from
#      pyclical/demos/plotting/plotting-env.yml.
#   2. Activates the environment.
#   3. On systems with a native GPU driver (detected via /dev/dri/card0),
#      removes the conda-forge mesalib package to prevent it from conflicting
#      with the system OpenGL driver.
#
# Do not run this script on ARM aarch64 (Asahi Linux): conda-forge VTK/Mayavi
# binaries are 4KB page-aligned and will segfault on Asahi's 16KB page system.
# See INSTALL.md for the ARM-specific venv procedure.

if [ "${BASH_SOURCE-}" = "$0" ]; then
    echo "Please run this script from Bash as follows:" >&2
    echo "  source $0" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.."; pwd)"

if ! command -v mamba >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
    echo "Error: neither mamba nor conda found on PATH." >&2
    echo "Install Miniforge (https://github.com/conda-forge/miniforge) and try again." >&2
    return 1
fi

echo "Creating/updating Conda environment from pyclical/demos/plotting/plotting-env.yml ..."
if command -v mamba >/dev/null 2>&1; then
    mamba env create -f "${REPO_ROOT}/pyclical/demos/plotting/plotting-env.yml" 2>/dev/null \
        || mamba env update -f "${REPO_ROOT}/pyclical/demos/plotting/plotting-env.yml"
else
    conda env create -f "${REPO_ROOT}/pyclical/demos/plotting/plotting-env.yml" 2>/dev/null \
        || conda env update -f "${REPO_ROOT}/pyclical/demos/plotting/plotting-env.yml"
fi

if command -v mamba >/dev/null 2>&1; then
    mamba activate pyclical-plotting
else
    conda activate pyclical-plotting
fi

# Remove conda-forge's Mesa library if a native GPU driver is present.
# The presence of /dev/dri/card0 indicates a real GPU with its own OpenGL
# driver. The conda-forge mesalib conflicts with that driver and must be
# removed to prevent rendering failures.
if [ -e /dev/dri/card0 ]; then
    echo "Native GPU detected (/dev/dri/card0 exists)."
    if conda list -n pyclical-plotting "^mesalib$" | grep -q mesalib; then
        echo "Removing conda-forge mesalib to avoid OpenGL conflict with system driver."
        conda remove -n pyclical-plotting --force mesalib -y
    else
        echo "conda-forge mesalib is not present in the environment; no removal needed."
    fi
else
    echo "No native GPU device detected; keeping conda-forge mesalib for software rendering."
fi

echo ""
echo "Environment 'pyclical-plotting' is ready."
echo "Next steps (from the repository root, with this environment active):"
echo "  make -f admin/Makefile.common bootstrap  # git clone only, not needed for tarballs"
echo "  make -C pyclical -j\$(($(nproc)/2))"
echo "  source pyclical/demos/plotting/export-plotting-vars.sh"
