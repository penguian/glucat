AUTOMAKE_OPTIONS = foreign

if make_pyclical
plotdir = $(demodir)/plotting
dist_plot_DATA = setup-plotting-env.sh export-plotting-vars.sh plotting-env.yml \
    plotting_demo_mayavi.py plotting_demo_dialog.py
endif

CLEANFILES = *.gcno *.gcda
