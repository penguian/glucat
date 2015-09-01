####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg06.cpp peg06.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test06

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test06_LDFLAGS = $(LDFLAGS)

test06_SOURCES = peg06.cpp
test06_LDADD = $(LDADD)
test06_DEPENDENCIES = ../test/driver.h peg06.h
