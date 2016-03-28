####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg00.cpp peg00.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test00
check_DATA = test00.out
test00.out: test00
	../test/test.sh 00 > test00.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test00_LDFLAGS = $(LDFLAGS)

test00_SOURCES = peg00.cpp
test00_LDADD = $(LDADD)
test00_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg00.h

clean:
	test -z peg00.o    || rm -f peg00.o
	test -z test00     || rm -f test00
	test -z test00.out || rm -f test00.out
