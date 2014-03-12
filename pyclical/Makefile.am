EXTRA_DIST = demos/pyclical_demo.py demos/pyclical_demo.out Makefile.am

AUTOMAKE_OPTIONS = foreign

ext_name = PyClical

# Setup and clean options depend on availability of Cython
if use_cython
  SETUP_SOURCE_LIST = $(ext_name).pyx $(ext_name).pxd glucat.pxd $(ext_name).h
  SETUP_SOURCE_PARAM = source_pyx="$(ext_name).pyx"
  SETUP = setup
  CLEAN_CYTHON = $(ext_name).cpp $(ext_name).html
else
  SETUP_SOURCE_LIST = $(ext_name).cpp $(ext_name).h
  SETUP_SOURCE_PARAM = source_cpp="$(ext_name).cpp"
  SETUP = setup_nocython
  CLEAN_CYTHON =
endif

# C++ compiler flags
glucat_extra_cxxflags_pyclical = @glucat_extra_cxxflags_pyclical@
EXTCXXFLAGS = $(glucat_extra_cxxflags_pyclical) $(CXXFLAGS)
# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)

pyexec_DATA = $(ext_name).so

all: $(pyexec_DATA)

$(ext_name).so: $(SETUP_SOURCE_LIST) $(SETUP).py setup_ext.py ../glucat/glucat.h ../glucat/glucat_imp.h
	ext_name="$(ext_name)" $(SETUP_SOURCE_PARAM) CXXFLAGS="$(EXTCXXFLAGS)" INCLUDES="$(INCLUDES)" LDFLAGS="$(LDFLAGS)" $(PYTHON) $(SETUP).py build_ext --inplace

clean:
	rm -f $(ext_name).so $(CLEAN_CYTHON)
	rm -rf build
