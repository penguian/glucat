EXTRA_DIST = pyclical_demo.py pyclical_demo.out Makefile.am

AUTOMAKE_OPTIONS = foreign

# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)

all: PyCliCal.so PyCliCal.html

PyCliCal.so PyCliCal.html: setup.py PyCliCal.pyx PyCliCal.h
	rm -f PyCliCal.so PyCliCal.cpp PyCliCal.html
	CXXFLAGS="$(CXXFLAGS)" INCLUDES="$(INCLUDES)" LDFLAGS="$(LDFLAGS)" python setup.py build_ext --inplace

clean:
	rm -f PyCliCal.so PyCliCal.cpp PyCliCal.html
	rm -rf build