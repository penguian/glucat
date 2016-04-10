EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test03

check: test03.out
test03.out: test03
	../test/test.sh 03 > test03.out

noinst_HEADERS = peg03.h
test03_SOURCES = peg03.cpp
test03_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test03_LDFLAGS = $(LDFLAGS)

CLEANFILES = test03.out
