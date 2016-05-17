AUTOMAKE_OPTIONS = foreign

if make_pyclical

DEMOS =                      pyclical_demo.py                     pyclical_tutorial_0_4_transcendental.py  pyclical_tutorial_1_5_conformal.py \
build_pyclical_notebooks.py  pyclical_notebook_utils.py           pyclical_tutorial_1_0_plane.py           pyclical_tutorial_2_0_wedge_product.py \
clifford_demo.py             pyclical_tutorial_0_0_notation.py    pyclical_tutorial_1_1_complex.py         pyclical_tutorials.py \
plotting_demo_dialog.py      pyclical_tutorial_0_1_index_sets.py  pyclical_tutorial_1_2_space.py           pyclical_tutorial_utils.py \
plotting_demo_mayavi.py      pyclical_tutorial_0_2_operations.py  pyclical_tutorial_1_3_lorentz.py         sqrt_log_demo.py \
plotting_demo.py             pyclical_tutorial_0_3_functions.py   pyclical_tutorial_1_4_fourth.py 

DEMOSOUT = $(DEMOS) pyclical_demo.out  sqrt_log_demo.out

dist_demo_DATA = $(DEMOSOUT)

NOTEBOOKS = clifford_demo.ipynb            pyclical_demo.ipynb                    sqrt_log_demo.ipynb \
pyclical_tutorial_0_0_notation.ipynb       pyclical_tutorial_0_1_index_sets.ipynb pyclical_tutorial_0_2_operations.ipynb   pyclical_tutorial_0_3_functions.ipynb \
pyclical_tutorial_0_4_transcendental.ipynb pyclical_tutorial_1_0_plane.ipynb      pyclical_tutorial_1_1_complex.ipynb      pyclical_tutorial_1_2_space.ipynb \
pyclical_tutorial_1_3_lorentz.ipynb        pyclical_tutorial_1_4_fourth.ipynb     pyclical_tutorial_1_5_conformal.ipynb    pyclical_tutorial_2_0_wedge_product.ipynb

nodist_demo_DATA = $(NOTEBOOKS)

$(NOTEBOOKS):
	PYTHONPATH=..:$(PYTHONPATH) python build_pyclical_notebooks.py

clean-local:
	rm -f *.pyc *.ipynb

endif
