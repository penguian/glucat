# Use conda and pip to install PyQt5 and Mayavi

# Run from bash using source
if [ "${BASH_SOURCE-}" = "$0" ]; then
  echo "Please run this script from Bash as follows: \$ source $0" >&2
  exit 1
fi
conda create -y -n mayavi
conda activate mayavi
conda install 'python=3.9' -y -c conda-forge
pip uninstall PyQt5
pip install PyQt5
conda install mayavi -y -c conda-forge
conda uninstall mesalib -y
conda install cython -y -c conda-forge
conda install matplotlib -y -c conda-forge
conda install qd -y -c conda-forge
conda uninstall mayavi -y
conda install mayavi -y -c conda-forge
