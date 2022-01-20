EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test01

check: test01.out
test01.out: test01
	../test/test.sh 01 > test01.out

noinst_HEADERS = peg01.h
test01_SOURCES = peg01.cpp
test01_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test01_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test01.out
