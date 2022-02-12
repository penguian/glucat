EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test10

check: test10.out
test10.out: test10
	../test/test.sh 10 > test10.out

noinst_HEADERS = peg10.h
test10_SOURCES = peg10.cpp
test10_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test10_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test10.out
