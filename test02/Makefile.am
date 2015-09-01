####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg02.cpp peg02.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test02

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test02_LDFLAGS = $(LDFLAGS)

test02_SOURCES = peg02.cpp
test02_LDADD = $(LDADD)
test02_DEPENDENCIES = ../test/driver.h peg02.h
