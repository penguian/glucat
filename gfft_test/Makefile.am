####### kdevelop will overwrite this part!!! (begin)##########
INCLUDES = $(all_includes)
EXTRA_DIST = Makefile.am gfft_test.cpp gfft_test.h
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
noinst_PROGRAMS = gfft_test

# set the include path found by configure
gfft_test_INCLUDES = $(all_includes)
# the library search path.
LDFLAGS =
gfft_test_LDFLAGS = $(all_libraries)

gfft_test_SOURCES = gfft_test.cpp
gfft_test_LDADD = $(LDADD)
gfft_test_DEPENDENCIES = gfft_test.h
