EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = squaring

noinst_HEADERS = squaring.h

squaring_SOURCES = squaring.cpp
squaring_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
squaring_LDFLAGS = $(USER_LDFLAGS)
