####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg10.cpp peg10.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test10
check_DATA = test10.out
test10.out: test10
	../test/test.sh 10 > test10.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test10_LDFLAGS = $(LDFLAGS)

test10_SOURCES = peg10.cpp
test10_LDADD = $(LDADD)
test10_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg10.h

clean:
	test -z peg10.o    || rm -f peg10.o
	test -z test10     || rm -f test10
	test -z test10.out || rm -f test10.out
