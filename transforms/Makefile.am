####### kdevelop will overwrite this part!!! (begin)##########
INCLUDES = $(all_includes)
EXTRA_DIST = Makefile.am transforms.h transforms.h 
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = transforms

# set the include path found by configure
transforms_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
transforms_LDFLAGS = $(all_libraries)

transforms_SOURCES = transforms.cpp
transforms_LDADD = $(LDADD)
transforms_DEPENDENCIES = driver.h transforms.h
