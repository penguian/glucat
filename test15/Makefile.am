####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg15.cpp peg15.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test15

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test15_LDFLAGS = $(LDFLAGS)

test15_SOURCES = peg15.cpp
test15_LDADD = $(LDADD)
test15_DEPENDENCIES = ../test/driver.h peg15.h
