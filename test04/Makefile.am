EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test04

check: test04.out
test04.out: test04
	../test/test.sh 04 > test04.out

noinst_HEADERS = peg04.h
test04_SOURCES = peg04.cpp
test04_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test04_LDFLAGS = $(LDFLAGS)

CLEANFILES = test04.out
