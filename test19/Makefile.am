EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test19

check: test19.out
test19.out: test19
	../test/test.sh 19 > test19.out

noinst_HEADERS = peg19.h
test19_SOURCES = test19.cpp
test19_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
test19_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test19.out *.gcno *.gcda
