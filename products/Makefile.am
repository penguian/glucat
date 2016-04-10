EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = products

noinst_HEADERS = products.h

products_SOURCES = products.cpp
products_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
products_LDFLAGS = $(LDFLAGS)
