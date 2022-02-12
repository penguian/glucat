EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test00

check: test00.out
test00.out: test00
	../test/test.sh 00 > test00.out

noinst_HEADERS = peg00.h
test00_SOURCES = peg00.cpp
test00_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test00_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test00.out
