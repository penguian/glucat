####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg07.cpp peg07.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test07
check_DATA = test07.out
test07.out: test07
	../test/test.sh 07 > test07.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test07_LDFLAGS = $(LDFLAGS)

test07_SOURCES = peg07.cpp
test07_LDADD = $(LDADD)
test07_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg07.h

clean:
	rm -f peg07.o test07 test07.out
