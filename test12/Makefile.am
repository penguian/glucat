####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = ../test/driver.h peg12.cpp peg12.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test12

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test12_LDFLAGS = $(LDFLAGS)

test12_SOURCES = peg12.cpp
test12_LDADD = $(LDADD)
test12_DEPENDENCIES = ../test/driver.h peg12.h
