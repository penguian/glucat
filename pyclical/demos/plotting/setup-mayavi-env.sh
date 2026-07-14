#!/usr/bin/env bash
# Set up the pyclical-mayavi Conda environment for the Mayavi plotting demo.
#
# Run from the glucat repository root with:
#   source pyclical/demos/plotting/setup-mayavi-env.sh
#
# This script:
#   1. Creates (or updates) the pyclical-mayavi Conda environment from
#      pyclical/demos/plotting/mayavi-env.yml.
#   2. Activates the environment.
#   3. On systems with a native GPU driver (detected via /dev/dri/card0),
#      removes the conda-forge mesalib package to prevent it from conflicting
#      with the system OpenGL driver.
#
# Do NOT run this script on ARM aarch64 (Asahi Linux): conda-forge VTK/Mayavi
# binaries are 4KB page-aligned and will segfault on Asahi's 16KB page kernel.
# Use setup-pyvista-env.sh instead.

if [ "${BASH_SOURCE-}" = "$0" ]; then
    echo "Please run this script from Bash as follows:" >&2
    echo "  source $0" >&2
    exit 1
fi

_ARCH=$(uname -m)
if [ "$_ARCH" = "aarch64" ]; then
    echo "Error: setup-mayavi-env.sh is for x86-64 systems only." >&2
    echo "conda-forge Mayavi/VTK binaries fail on ARM64 16KB page size kernels." >&2
    echo "Please use setup-pyvista-env.sh for ARM64 / Asahi Linux." >&2
    return 1 2>/dev/null || exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.."; pwd)"

if ! command -v mamba >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
    echo "Error: neither mamba nor conda found on PATH." >&2
    echo "Install Miniforge (https://github.com/conda-forge/miniforge) and try again." >&2
    return 1
fi

# Select mamba if available, otherwise fall back to conda.
if command -v mamba >/dev/null 2>&1; then
    _CONDACMD=mamba
else
    _CONDACMD=conda
fi

echo "Creating/updating Conda environment from pyclical/demos/plotting/mayavi-env.yml ..."
$_CONDACMD env create -f "${REPO_ROOT}/pyclical/demos/plotting/mayavi-env.yml" 2>/dev/null \
    || $_CONDACMD env update -f "${REPO_ROOT}/pyclical/demos/plotting/mayavi-env.yml"

$_CONDACMD activate pyclical-mayavi

# Remove conda-forge's Mesa library if a native GPU driver is present.
if [ -e /dev/dri/card0 ]; then
    echo "Native GPU detected (/dev/dri/card0 exists)."
    if $_CONDACMD list -n pyclical-mayavi "^mesalib$" | grep -q mesalib; then
        echo "Removing conda-forge mesalib to avoid OpenGL conflict with system driver."
        $_CONDACMD remove -n pyclical-mayavi --force mesalib -y
    else
        echo "conda-forge mesalib is not present in the environment; no removal needed."
    fi
else
    echo "No native GPU device detected; keeping conda-forge mesalib for software rendering."
fi

echo ""
echo "Environment 'pyclical-mayavi' is ready."
echo "Next steps (from the repository root, with this environment active):"
echo "  make -f admin/Makefile.common bootstrap  # git clone only, not needed for tarballs"
echo "  ./configure"
echo "  make -C pyclical -j\$(( \$(nproc) > 1 ? \$(nproc)/2 : 1 ))"
echo "  source pyclical/demos/plotting/export-mayavi-vars.sh"
