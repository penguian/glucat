EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test11

check: test11.out
test11.out: test11
	../test/test.sh 11 > test11.out

noinst_HEADERS = peg11.h
test11_SOURCES = peg11.cpp
test11_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test11_LDFLAGS = $(LDFLAGS)

CLEANFILES = test11.out
