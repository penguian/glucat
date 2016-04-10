EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test07

check: test07.out
test07.out: test07
	../test/test.sh 07 > test07.out

noinst_HEADERS = peg07.h
test07_SOURCES = peg07.cpp
test07_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test07_LDFLAGS = $(LDFLAGS)

CLEANFILES = test07.out
