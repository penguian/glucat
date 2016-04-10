EXTRA_DIST = Makefile.am eg8.txt

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test08

check_DATA = eg8.txt

check: test08.out
test08.out: test08 eg8.txt
	../test/test.sh 08 > test08.out

noinst_HEADERS = peg08.h
test08_SOURCES = peg08.cpp
test08_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test08_LDFLAGS = $(LDFLAGS)

CLEANFILES = test08.out
