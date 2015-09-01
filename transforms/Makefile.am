####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = Makefile.am transforms.h transforms.h 
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = transforms

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
transforms_LDFLAGS = $(LDFLAGS)

transforms_SOURCES = transforms.cpp
transforms_LDADD = $(LDADD)
transforms_DEPENDENCIES = transforms.h
