####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg02.cpp peg02.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test02
check_DATA = test02.out
test02.out: test02
	../test/test.sh 02 > test02.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test02_LDFLAGS = $(LDFLAGS)

test02_SOURCES = peg02.cpp
test02_LDADD = $(LDADD)
test02_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg02.h

clean:
	rm -f peg02.o test02 test02.out
