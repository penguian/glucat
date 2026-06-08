EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = expressions

noinst_HEADERS = expressions.h

expressions_SOURCES = expressions.cpp
expressions_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
expressions_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = *.gcno *.gcda
