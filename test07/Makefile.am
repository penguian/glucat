####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg07.cpp peg07.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test07

# set the include path found by configure
test07_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
test07_LDFLAGS = $(all_libraries)

test07_SOURCES = peg07.cpp
test07_LDADD = $(LDADD)
test07_DEPENDENCIES = ../test/driver.h peg07.h
