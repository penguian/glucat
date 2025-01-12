# Run using source
# Use conda and pip to install PyQt5 and Mayavi

conda create -y -n mayavi
conda activate mayavi
conda install 'python=3.9' -y -c conda-forge
pip install PyQt5
conda install mayavi -y -c conda-forge
conda uninstall mesalib -y
