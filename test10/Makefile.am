####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg10.cpp peg10.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test10

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test10_LDFLAGS = $(LDFLAGS)

test10_SOURCES = peg10.cpp
test10_LDADD = $(LDADD)
test10_DEPENDENCIES = ../test/driver.h peg10.h
