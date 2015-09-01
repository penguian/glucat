####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg14.cpp peg14.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test14

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test14_LDFLAGS = $(LDFLAGS)

test14_SOURCES = peg14.cpp
test14_LDADD = $(LDADD)
test14_DEPENDENCIES = ../test/driver.h peg14.h
