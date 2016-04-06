####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = eg8.txt peg08.cpp peg08.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test08
check_DATA = test08.out eg8.txt
test08.out: test08 eg8.txt
	../test/test.sh 08 > test08.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test08_LDFLAGS = $(LDFLAGS)

test08_SOURCES = peg08.cpp
test08_LDADD = $(LDADD)
test08_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg08.h

clean:
	rm -f peg08.o test08 test08.out
