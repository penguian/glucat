####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg00.cpp peg00.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test00

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test00_LDFLAGS = $(LDFLAGS)

test00_SOURCES = peg00.cpp
test00_LDADD = $(LDADD)
test00_DEPENDENCIES = ../test/driver.h peg00.h
