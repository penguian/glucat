####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = Makefile.am gfft_test.cpp gfft_test.h
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = gfft_test

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
gfft_test_LDFLAGS = $(LDFLAGS)

gfft_test_SOURCES = gfft_test.cpp
gfft_test_LDADD = $(LDADD)
gfft_test_DEPENDENCIES = gfft_test.h
