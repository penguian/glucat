####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg05.cpp peg05.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test05

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test05_LDFLAGS = $(LDFLAGS)

test05_SOURCES = peg05.cpp
test05_LDADD = $(LDADD)
test05_DEPENDENCIES = ../test/driver.h peg05.h
