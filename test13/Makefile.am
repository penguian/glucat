EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test13

check: test13.out
test13.out: test13
	../test/test.sh 13 > test13.out

noinst_HEADERS = peg13.h
test13_SOURCES = peg13.cpp
test13_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test13_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test13.out
