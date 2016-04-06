####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg13.cpp peg13.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test13
check_DATA = test13.out
test13.out: test13
	../test/test.sh 13 > test13.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test13_LDFLAGS = $(LDFLAGS)

test13_SOURCES = peg13.cpp
test13_LDADD = $(LDADD)
test13_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg13.h

clean:
	rm -f peg13.o test13 test13.out
