####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg08.cpp peg08.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test08

# set the include path found by configure
test08_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
test08_LDFLAGS = $(all_libraries)

test08_SOURCES = peg08.cpp
test08_LDADD = $(LDADD)
test08_DEPENDENCIES = ../test/driver.h peg08.h
