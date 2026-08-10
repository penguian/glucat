#!/usr/bin/env bash
# Set up the unified Python environment for PyClical (testing, linting, notebooks, and PyVista plotting).
#
# Run from the glucat repository root with:
#   source pyclical/setup-env.sh
#

if [ "${BASH_SOURCE-}" = "$0" ]; then
    echo "Please run this script from Bash as follows:" >&2
    echo "  source $0" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.."; pwd)"
_ARCH=$(uname -m)

if [ "$_ARCH" = "aarch64" ]; then
    echo "Configuring PyClical environment for ARM64 (aarch64) system site-packages ..."
    if ! python3 -c 'import vtkmodules' 2>/dev/null; then
        echo "Error: python3-vtk is not installed in system site-packages." >&2
        echo "The 16KB page-aligned system VTK RPM is required on ARM64 / Asahi Linux." >&2
        echo "Please run: sudo dnf install python3-vtk" >&2
        return 1 2>/dev/null || exit 1
    fi

    mkdir -p "${REPO_ROOT}/.venvs"
    python3 -m venv --system-site-packages "${REPO_ROOT}/.venvs/pyclical-env"
    source "${REPO_ROOT}/.venvs/pyclical-env/bin/activate"

    echo "Installing PyVista without bundled VTK (reusing system python3-vtk RPM)..."
    pip install --no-deps pyvista

    grep -v '^pyvista$' "${REPO_ROOT}/pyclical/requirements.txt" > "${REPO_ROOT}/.venvs/arm_requirements.txt"
    pip install -r "${REPO_ROOT}/.venvs/arm_requirements.txt"
    rm -f "${REPO_ROOT}/.venvs/arm_requirements.txt"

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
    echo "Virtual environment '.venvs/pyclical-env' is ready."
else
    if command -v mamba >/dev/null 2>&1; then
        _CONDACMD=mamba
    elif command -v conda >/dev/null 2>&1; then
        _CONDACMD=conda
    else
        _CONDACMD=""
    fi

    if [ -n "$_CONDACMD" ]; then
        echo "Creating/updating Conda environment from pyclical/environment.yml ..."
        $_CONDACMD env update -f "${REPO_ROOT}/pyclical/environment.yml" 2>/dev/null \
            || $_CONDACMD env create -f "${REPO_ROOT}/pyclical/environment.yml"

        if [ "$_CONDACMD" = "mamba" ]; then
            eval "$(mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
        else
            eval "$($_CONDACMD shell.bash hook 2>/dev/null)"
        fi
        $_CONDACMD activate pyclical-env
        echo ""
        echo "Conda environment 'pyclical-env' is ready."
    else
        echo "Creating standard virtual environment at .venvs/pyclical-env ..."
        mkdir -p "${REPO_ROOT}/.venvs"
        python3 -m venv "${REPO_ROOT}/.venvs/pyclical-env"
        source "${REPO_ROOT}/.venvs/pyclical-env/bin/activate"
        pip install -r "${REPO_ROOT}/pyclical/requirements.txt"
        echo ""
        echo "Virtual environment '.venvs/pyclical-env' is ready."
    fi
fi

# Export environment variables for PyVista plotting and PyClical imports
export PYTHONPATH="${REPO_ROOT}/pyclical:${REPO_ROOT}/pyclical/demos/plotting:${PYTHONPATH}"
export QT_API="pyside6"
export QT_QPA_PLATFORM="xcb"

# Register Jupyter Kernel if PyClical build targets exist
if [ -f "${REPO_ROOT}/pyclical/Makefile" ]; then
    make -C "${REPO_ROOT}/pyclical" install-pyclical-kernel 2>/dev/null || true
fi

echo "Next steps (from the repository root, with this environment active):"
echo "  ./configure PYTHON=\$(which python3)"
echo "  make -C pyclical -j\$(( \$(nproc) > 1 ? \$(nproc)/2 : 1 ))"
echo "  ./verify_all.py --fast --coverage --python"
