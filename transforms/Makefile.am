####### kdevelop will overwrite this part!!! (begin)##########
INCLUDES = $(all_includes)
EXTRA_DIST = Makefile.am driver.cpp driver.h transforms.h 
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = glucat

# set the include path found by configure
glucat_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
glucat_LDFLAGS = $(all_libraries)

glucat_SOURCES = driver.cpp
glucat_LDADD = $(LDADD)
glucat_DEPENDENCIES = driver.h transforms.h
