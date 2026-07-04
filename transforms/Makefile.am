EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = transforms

noinst_HEADERS = transforms.h

transforms_SOURCES = transforms.cpp
transforms_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes) $(GLUCAT_CPPFLAGS)
# the library search path.
transforms_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = *.gcno *.gcda
