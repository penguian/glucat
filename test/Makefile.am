####### kdevelop will overwrite this part!!! (begin)##########

INCLUDES = $(all_includes)

glucat_SOURCES = driver.cpp 

EXTRA_DIST = driver.cpp driver.h eg8.txt peg01.h peg02.h peg03.h peg04.h peg05.h peg06.h peg07.h peg08.h peg09.h peg10.h peg11.h peg12.h peg13.h peg14.h peg15.h peg16.h squaring.h eg3.res test-6.out test.out 

####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = glucat
AM_CXXFLAGS = -ftemplate-depth-32 -fexceptions -ansi -Wold-style-cast -Woverloaded-virtual -Wall-Wno-return-type

# Uncomment for normal test build - excluding Compaq Alpha:
CXXFLAGS = -O1 -finline-limit=400
# Uncomment for normal test build on Compaq Alpha:
# CXXFLAGS = -O1 -finline-limit=400 -mieee
# Uncomment for debugging - excluding Compaq Alpha:
# CXXFLAGS = -O1 -finline-limit=400 -g3
# Uncomment for speed test build - excluding Compaq Alpha:
# CXXFLAGS = -O2 -finline-limit=800 -finline-functions
# Uncomment for speed test build on Compaq Alpha:
# CXXFLAGS = -O2 -finline-limit=800 -finline-functions -mieee

# set the include path found by configure
glucat_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
glucat_LDFLAGS = $(all_libraries)
