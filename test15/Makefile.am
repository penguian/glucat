####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg15.cpp peg15.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test15
check_DATA = test15.out
test15.out: test15
	../test/test.sh 15 > test15.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test15_LDFLAGS = $(LDFLAGS)

test15_SOURCES = peg15.cpp
test15_LDADD = $(LDADD)
test15_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg15.h

clean:
	test -z peg15.o    || rm -f peg15.o
	test -z test15     || rm -f test15
	test -z test15.out || rm -f test15.out
