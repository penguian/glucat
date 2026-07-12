# Define environment variables for Mayavi plotting demos

# Run from bash using source
if [ "${BASH_SOURCE-}" = "$0" ]; then
  echo "Please run this script from Bash as follows: \$ source $0" >&2
  exit 1
fi
export PYTHONPATH="$(cd ..;pwd):$PYTHONPATH"
export QT_API="pyqt5"
export QT_QPA_PLATFORM="xcb"
