EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test02

check: test02.out
test02.out: test02
	../test/test.sh 02 > test02.out

noinst_HEADERS = peg02.h
test02_SOURCES = peg02.cpp
test02_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test02_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test02.out
