AUTOMAKE_OPTIONS = foreign

SUBDIRS = plotting

if make_pyclical

DEMOS =  build_pyclical_notebooks.py    clifford_demo.py                         m_theory_demo.py \
pyclical_demo.py                        pyclical_notebook_utils.py               sqrt_log_demo.py \
pyclical_tutorial_0_0_notation.py       pyclical_tutorial_0_1_index_sets.py      pyclical_tutorial_0_2_operations.py \
pyclical_tutorial_0_3_functions.py      pyclical_tutorial_0_4_transcendental.py  pyclical_tutorial_1_0_plane.py \
pyclical_tutorial_1_1_complex.py        pyclical_tutorial_1_2_space.py           pyclical_tutorial_1_3_lorentz.py \
pyclical_tutorial_1_4_fourth.py         pyclical_tutorial_1_5_conformal.py       pyclical_tutorial_2_0_wedge_product.py \
pyclical_tutorial_utils.py              pyclical_tutorials.py                    versor_demo.py

DEMOSOUT = $(DEMOS) clifford_demo.out pyclical_demo.out sqrt_log_demo.out

dist_demo_DATA = $(DEMOSOUT)

NOTEBOOKS = clifford_demo.ipynb        m_theory_demo.ipynb                         pyclical_demo.ipynb \
sqrt_log_demo.ipynb                     versor_demo.ipynb \
pyclical_tutorial_0_0_notation.ipynb   pyclical_tutorial_0_1_index_sets.ipynb      pyclical_tutorial_0_2_operations.ipynb \
pyclical_tutorial_0_3_functions.ipynb  pyclical_tutorial_0_4_transcendental.ipynb  pyclical_tutorial_1_0_plane.ipynb \
pyclical_tutorial_1_1_complex.ipynb    pyclical_tutorial_1_2_space.ipynb           pyclical_tutorial_1_3_lorentz.ipynb \
pyclical_tutorial_1_4_fourth.ipynb     pyclical_tutorial_1_5_conformal.ipynb       pyclical_tutorial_2_0_wedge_product.ipynb

nodist_demo_DATA = $(NOTEBOOKS)

notebooks.stamp: build_pyclical_notebooks.py
	PYTHONPATH="..:$(PYTHONPATH)" $(PYTHON) build_pyclical_notebooks.py
	touch notebooks.stamp

$(NOTEBOOKS): notebooks.stamp
	@true

clean-local:
	rm -rf *.pyc *.ipynb __pycache__ notebooks.stamp

endif

CLEANFILES = *.gcno *.gcda
