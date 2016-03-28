####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg06.cpp peg06.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test06
check_DATA = test06.out
test06.out: test06
	../test/test.sh 06 > test06.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test06_LDFLAGS = $(LDFLAGS)

test06_SOURCES = peg06.cpp
test06_LDADD = $(LDADD)
test06_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg06.h

clean:
	test -z peg06.o    || rm -f peg06.o
	test -z test06     || rm -f test06
	test -z test06.out || rm -f test06.out
