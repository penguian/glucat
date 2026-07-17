# Define environment variables for PyVista plotting demos

# Run from the glucat repository root using source:
#   source pyclical/demos/plotting/export-pyvista-vars.sh
if [ "${BASH_SOURCE-}" = "$0" ]; then
  echo "Please run this script from Bash as follows: \$ source $0" >&2
  exit 1
fi
export PYTHONPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/../..";pwd):$PYTHONPATH"
export QT_API="pyside6"
export QT_QPA_PLATFORM="xcb"
