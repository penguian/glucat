EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test18

check: test18.out
test18.out: test18
	../test/test.sh 18 > test18.out

noinst_HEADERS = peg18.h
test18_SOURCES = peg18.cpp
test18_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
test18_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test18.out *.gcno *.gcda
