####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg01.cpp peg01.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test01

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test01_LDFLAGS = $(LDFLAGS)

test01_SOURCES = peg01.cpp
test01_LDADD = $(LDADD)
test01_DEPENDENCIES = ../test/driver.h peg01.h
