####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = squaring.cpp squaring.h tuning.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = squaring

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
squaring_LDFLAGS = $(LDFLAGS)

squaring_SOURCES = squaring.cpp
squaring_LDADD = $(LDADD)
squaring_DEPENDENCIES = squaring.h
