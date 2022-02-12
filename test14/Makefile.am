EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test14

check: test14.out
test14.out: test14
	../test/test.sh 14 > test14.out

noinst_HEADERS = peg14.h
test14_SOURCES = peg14.cpp
test14_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test14_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test14.out
