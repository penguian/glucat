####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = products.cpp products.h tuning.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = products

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
products_LDFLAGS = $(LDFLAGS)

products_SOURCES = products.cpp
products_LDADD = $(LDADD)
products_DEPENDENCIES = products.h
