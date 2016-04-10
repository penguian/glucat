EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = gfft_test

noinst_HEADERS = gfft_test.h

gfft_test_SOURCES = gfft_test.cpp
gfft_test_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
gfft_test_LDFLAGS = $(LDFLAGS)
