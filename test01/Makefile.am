####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = driver.h peg01.cpp peg01.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = test01

# set the include path found by configure
test01_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
test01_LDFLAGS = $(all_libraries)

test01_SOURCES = peg01.cpp
test01_LDADD = $(LDADD)
test01_DEPENDENCIES = ../test/driver.h peg01.h
