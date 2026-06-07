EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = gfft_test

noinst_HEADERS = gfft_test.h

gfft_test_SOURCES = gfft_test.cpp
gfft_test_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
gfft_test_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = *.gcno *.gcda
