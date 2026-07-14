AUTOMAKE_OPTIONS = foreign

if make_pyclical
plotdir = $(demodir)/plotting
dist_plot_DATA = \
    setup-mayavi-env.sh \
    setup-pyvista-env.sh \
    export-mayavi-vars.sh \
    export-pyvista-vars.sh \
    mayavi-env.yml \
    pyvista-env.yml \
    requirements-pyvista.txt
dist_plot_SCRIPTS = \
    plotting_demo_pyvista.py \
    plotting_demo_pyvista_dialog.py
endif

CLEANFILES = *.gcno *.gcda
