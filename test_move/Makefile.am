EXTRA_DIST = Makefile.am

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = test_move

check: test_move.out
test_move.out: test_move
	./test_move > test_move.out

test_move_SOURCES = test_move.cpp
test_move_LDADD = $(LDADD)

# set the include path found by configure
AM_CPPFLAGS = $(all_includes)
# the library search path.
test_move_LDFLAGS = $(USER_LDFLAGS)

CLEANFILES = test_move.out
