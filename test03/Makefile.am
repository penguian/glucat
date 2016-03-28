####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg03.cpp peg03.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test03
check_DATA = test03.out
test03.out: test03
	../test/test.sh 03 > test03.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test03_LDFLAGS = $(LDFLAGS)

test03_SOURCES = peg03.cpp
test03_LDADD = $(LDADD)
test03_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg03.h

clean:
	test -z peg03.o    || rm -f peg03.o
	test -z test03     || rm -f test03
	test -z test03.out || rm -f test03.out
