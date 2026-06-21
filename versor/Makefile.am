EXTRA_DIST = Makefile.am 

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = versor

noinst_HEADERS = versor.h

versor_SOURCES = versor.cpp
versor_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
versor_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = *.gcno *.gcda
