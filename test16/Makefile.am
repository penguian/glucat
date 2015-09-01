####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg16.cpp peg16.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test16

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test16_LDFLAGS = $(LDFLAGS)

test16_SOURCES = peg16.cpp
test16_LDADD = $(LDADD)
test16_DEPENDENCIES = ../test/driver.h peg16.h
