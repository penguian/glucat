####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg14.cpp peg14.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test14
check_DATA = test14.out
test14.out: test14
	../test/test.sh 14 > test14.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test14_LDFLAGS = $(LDFLAGS)

test14_SOURCES = peg14.cpp
test14_LDADD = $(LDADD)
test14_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg14.h

clean:
	rm -f peg14.o test14 test14.out
