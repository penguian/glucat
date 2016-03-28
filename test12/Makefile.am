####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg12.cpp peg12.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test12
check_DATA = test12.out
test12.out: test12
	../test/test.sh 12 > test12.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test12_LDFLAGS = $(LDFLAGS)

test12_SOURCES = peg12.cpp
test12_LDADD = $(LDADD)
test12_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg12.h

clean:
	test -z peg12.o    || rm -f peg12.o
	test -z test12     || rm -f test12
	test -z test12.out || rm -f test12.out
