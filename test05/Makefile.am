EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test05

check: test05.out
test05.out: test05
	../test/test.sh 05 > test05.out

noinst_HEADERS = peg05.h
test05_SOURCES = peg05.cpp
test05_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test05_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test05.out
