AUTOMAKE_OPTIONS = foreign

if make_pyclical

ext_name = PyClical

EXTNOCYTHONCPP = $(ext_name)_nocython.cpp

EXTRA_DIST = demos/*.py demos/*.out Makefile.am $(EXTNOCYTHONCPP) test.out

# Setup options depend on availability of Cython
if use_cython
  SETUP_SOURCE_LIST = $(ext_name).pyx $(ext_name).pxd glucat.pxd $(ext_name).h
  SETUP_SOURCE_PARAM = source_pyx="$(ext_name).pyx"
  SETUP = setup
else
  SETUP_SOURCE_LIST = $(EXTNOCYTHONCPP) $(ext_name).h
  SETUP_SOURCE_PARAM = source_cpp="$(EXTNOCYTHONCPP)"
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

$(pyexec_DATA): $(SETUP_SOURCE_LIST) $(SETUP).py setup_ext.py
	ext_name="$(ext_name)" $(SETUP_SOURCE_PARAM) CXX="$(CXX)" CXXFLAGS="$(EXTCXXFLAGS)" AM_CPPFLAGS="$(EXTAM_CPPFLAGS)" LDFLAGS="$(EXTLDFLAGS)" $(PYTHON) $(SETUP).py build_ext --inplace

if use_cython
$(ext_name).cpp: $(pyexec_DATA)

$(EXTNOCYTHONCPP): $(ext_name).cpp
	cp $(ext_name).cpp $@
	$(MAKE) clean

endif

CLEAN_CYTHON = $(ext_name).cpp $(ext_name).html

CLEANFILES = $(pyexec_DATA) $(CLEAN_CYTHON)

clean-local:
	rm -f *.pyc
	rm -rf build

endif
