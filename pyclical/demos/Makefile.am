AUTOMAKE_OPTIONS = foreign

demodir = $(datarootdir)/pyclical/demos

NOTEBOOKS = clifford_demo.ipynb            pyclical_demo.ipynb                    sqrt_log_demo.ipynb \
pyclical_tutorial_0_0_notation.ipynb       pyclical_tutorial_0_1_index_sets.ipynb pyclical_tutorial_0_2_operations.ipynb   pyclical_tutorial_0_3_functions.ipynb \
pyclical_tutorial_0_4_transcendental.ipynb pyclical_tutorial_1_0_plane.ipynb      pyclical_tutorial_1_1_complex.ipynb      pyclical_tutorial_1_2_space.ipynb \
pyclical_tutorial_1_3_lorentz.ipynb        pyclical_tutorial_1_4_fourth.ipynb     pyclical_tutorial_1_5_conformal.ipynb    pyclical_tutorial_2_0_wedge_product.ipynb

demo_DATA = $(NOTEBOOKS) *.py

$(NOTEBOOKS):
	python build_pyclical_notebooks.py

clean-local:
	rm -f *.pyc *.ipynb
