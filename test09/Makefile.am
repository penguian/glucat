EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test09

check: test09.out
test09.out: test09
	../test/test.sh 09 > test09.out

noinst_HEADERS = peg09.h
test09_SOURCES = peg09.cpp
test09_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test09_LDFLAGS = $(LDFLAGS)

CLEANFILES = test09.out
