EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = products

noinst_HEADERS = products.h

products_SOURCES = products.cpp
products_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
products_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = *.gcno *.gcda
