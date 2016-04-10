EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test15

check: test15.out
test15.out: test15
	../test/test.sh 15 > test15.out

noinst_HEADERS = peg15.h
test15_SOURCES = peg15.cpp
test15_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test15_LDFLAGS = $(LDFLAGS)

CLEANFILES = test15.out
