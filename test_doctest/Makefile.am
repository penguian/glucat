EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test_doctest

check: test_doctest.out
test_doctest.out: test_doctest
	./test_doctest > test_doctest.out

test_doctest_SOURCES = main.cpp
test_doctest_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test_doctest_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test_doctest.out
