####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg09.cpp peg09.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test09

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test09_LDFLAGS = $(LDFLAGS)

test09_SOURCES = peg09.cpp
test09_LDADD = $(LDADD)
test09_DEPENDENCIES = ../test/driver.h peg09.h
