EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test17

check: test17.out
test17.out: test17
	../test/test.sh 17 > test17.out

noinst_HEADERS = peg17.h
test17_SOURCES = peg17.cpp
test17_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test17_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test17.out
