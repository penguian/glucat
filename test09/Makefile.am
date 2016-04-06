####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg09.cpp peg09.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test09
check_DATA = test09.out
test09.out: test09
	../test/test.sh 09 > test09.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test09_LDFLAGS = $(LDFLAGS)

test09_SOURCES = peg09.cpp
test09_LDADD = $(LDADD)
test09_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg09.h

clean:
	rm -f peg09.o test09 test09.out
