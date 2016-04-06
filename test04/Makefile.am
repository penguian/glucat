####### kdevelop will overwrite this part!!! (begin)##########
EXTRA_DIST = peg04.cpp peg04.h Makefile.am
####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
check_PROGRAMS = test04
check_DATA = test04.out
test04.out: test04
	../test/test.sh 04 > test04.out

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)
test04_LDFLAGS = $(LDFLAGS)

test04_SOURCES = peg04.cpp
test04_LDADD = $(LDADD)
test04_DEPENDENCIES = ../test/control.h ../test/driver.h ../test/try_catch.h ../test/tuning.h peg04.h

clean:
	rm -f peg04.o test04 test04.out
