####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg13.cpp peg13.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test13

# set the include path found by configure
test13_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
test13_LDFLAGS = $(all_libraries)

test13_SOURCES = peg13.cpp
test13_LDADD = $(LDADD)
test13_DEPENDENCIES = ../test/driver.h peg13.h
