EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test06

check: test06.out
test06.out: test06
	../test/test.sh 06 > test06.out

noinst_HEADERS = peg06.h
test06_SOURCES = peg06.cpp
test06_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test06_LDFLAGS = $(LDFLAGS)

CLEANFILES = test06.out
