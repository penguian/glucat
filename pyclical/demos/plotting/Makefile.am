AUTOMAKE_OPTIONS = foreign

if make_pyclical
plotdir = $(demodir)/plotting
dist_plot_DATA = \
    setup-mayavi-env.sh \
    export-mayavi-vars.sh \
    mayavi-env.yml
dist_plot_SCRIPTS = \
    plotting_demo_mayavi.py \
    plotting_demo_dialog.py \
    plotting_demo_pyvista.py \
    plotting_demo_pyvista_dialog.py
endif

CLEANFILES = *.gcno *.gcda *.pyc

clean-local:
	rm -rf __pycache__
