####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg05.cpp peg05.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test05
check_DATA = test05.out
test05.out: test05
	../test/test.sh 05 > test05.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test05_LDFLAGS = $(LDFLAGS)

test05_SOURCES = peg05.cpp
test05_LDADD = $(LDADD)
test05_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg05.h

clean:
	test -z peg05.o    || rm -f peg05.o
	test -z test05     || rm -f test05
	test -z test05.out || rm -f test05.out
