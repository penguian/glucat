EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test16

check: test16.out
test16.out: test16
	../test/test.sh 16 > test16.out

noinst_HEADERS = peg16.h
test16_SOURCES = peg16.cpp
test16_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test16_LDFLAGS = $(LDFLAGS)

CLEANFILES = test16.out
