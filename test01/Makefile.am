####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg01.cpp peg01.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test01
check_DATA = test01.out
test01.out: test01
	../test/test.sh 01 > test01.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test01_LDFLAGS = $(LDFLAGS)

test01_SOURCES = peg01.cpp
test01_LDADD = $(LDADD)
test01_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg01.h

clean:
	test -z peg01.o    || rm -f peg01.o
	test -z test01     || rm -f test01
	test -z test01.out || rm -f test01.out
