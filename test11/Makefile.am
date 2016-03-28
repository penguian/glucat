####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg11.cpp peg11.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test11
check_DATA = test11.out
test11.out: test11
	../test/test.sh 11 > test11.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test11_LDFLAGS = $(LDFLAGS)

test11_SOURCES = peg11.cpp
test11_LDADD = $(LDADD)
test11_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg11.h

clean:
	test -z peg11.o    || rm -f peg11.o
	test -z test11     || rm -f test11
	test -z test11.out || rm -f test11.out
