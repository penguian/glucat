AUTOMAKE_OPTIONS = foreign

if make_pyclical
ext_name = PyClical

EXTRA_DIST = demos/*.py demos/*.out Makefile.am $(ext_name)_nocython.cpp test.out

# Setup options depend on availability of Cython
if use_cython
  SETUP_SOURCE_LIST = $(ext_name).pyx $(ext_name).pxd glucat.pxd $(ext_name).h
  SETUP_SOURCE_PARAM = source_pyx="$(ext_name).pyx"
  SETUP = setup
else
  SETUP_SOURCE_LIST = $(ext_name)_nocython.cpp $(ext_name).h
  SETUP_SOURCE_PARAM = source_cpp="$(ext_name)_nocython.cpp"
  SETUP = setup_nocython
endif

# C++ compiler flags
EXTCXXFLAGS = $(glucat_extra_cxxflags_pyclical) $(CXXFLAGS)
# set the include path found by configure
EXTAM_CPPFLAGS = $(all_includes)
# the library search path.
EXTLDFLAGS = $(all_libraries)

pyexec_DATA = $(ext_name).so

all: $(pyexec_DATA)

$(ext_name).so: $(SETUP_SOURCE_LIST) $(SETUP).py setup_ext.py ../glucat/glucat.h ../glucat/glucat_imp.h
	ext_name="$(ext_name)" $(SETUP_SOURCE_PARAM) CXX="$(CXX)" CXXFLAGS="$(EXTCXXFLAGS)" AM_CPPFLAGS="$(EXTAM_CPPFLAGS)" LDFLAGS="$(EXTLDFLAGS)" $(PYTHON) $(SETUP).py build_ext --inplace

CLEAN_CYTHON = $(ext_name).cpp $(ext_name).html

CLEANFILES = $(ext_name).so $(CLEAN_CYTHON)

clean-local:
	rm -f *.pyc
	rm -rf build
endif
