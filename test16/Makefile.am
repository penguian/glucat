####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg16.cpp peg16.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test16
check_DATA = test16.out
test16.out: test16
	../test/test.sh 16 > test16.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test16_LDFLAGS = $(LDFLAGS)

test16_SOURCES = peg16.cpp
test16_LDADD = $(LDADD)
test16_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg16.h

clean:
	rm -f peg16.o test16 test16.out
