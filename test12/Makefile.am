EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test12

check: test12.out
test12.out: test12
	../test/test.sh 12 > test12.out

noinst_HEADERS = peg12.h
test12_SOURCES = peg12.cpp
test12_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test12_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test12.out
