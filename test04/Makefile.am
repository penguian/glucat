####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg04.cpp peg04.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test04

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test04_LDFLAGS = $(LDFLAGS)

test04_SOURCES = peg04.cpp
test04_LDADD = $(LDADD)
test04_DEPENDENCIES = ../test/driver.h peg04.h
