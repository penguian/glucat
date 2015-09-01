####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg03.cpp peg03.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test03

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test03_LDFLAGS = $(LDFLAGS)

test03_SOURCES = peg03.cpp
test03_LDADD = $(LDADD)
test03_DEPENDENCIES = ../test/driver.h peg03.h
