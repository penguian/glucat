# Use conda and pip to install PyQt5 and Mayavi

# Run from bash using source
if [ "${BASH_SOURCE-}" = "$0" ]; then
  echo "Please run this script from Bash as follows: \$ source $0" >&2
  exit 1
fi

conda create -n mayavi -y
conda activate mayavi
conda install 'python=3.9' -c conda-forge -y
conda install cython -c conda-forge -y
conda install matplotlib -c conda-forge -y
conda install qd -c conda-forge -y
conda install jupyter -c conda-forge -y
pip uninstall PyQt5
pip install PyQt5
conda install 'mayavi=4.8.2' -c conda-forge -y
conda install 'python=3.9' -c conda-forge -y
conda uninstall mesalib -y
