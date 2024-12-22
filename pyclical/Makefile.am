AUTOMAKE_OPTIONS = foreign

if make_pyclical

ext_name = PyClical

EXTCYTHON_SOURCE = $(ext_name).pyx $(ext_name).pxd glucat.pxd
EXTNOCYTHONCPP = $(ext_name)_nocython.cpp
EXTNOCYTHON_SOURCE = $(EXTNOCYTHONCPP)
EXT_HEADER = $(ext_name).h

EXTRA_DIST = Makefile.am $(EXTCYTHON_SOURCE) $(EXTNOCYTHON_SOURCE) $(EXT_HEADER) \
             setup_ext.py  setup_nocython.py  setup.py  test.py test.out

# Setup options depend on availability of Cython
if use_cython
  SETUP_SOURCE_LIST = $(EXTCYTHON_SOURCE) $(EXT_HEADER)
  SETUP_SOURCE_PARAM = source_pyx="$(ext_name).pyx"
  SETUP = setup
else
  SETUP_SOURCE_LIST = $(EXTNOCYTHON_SOURCE) $(EXT_HEADER)
  SETUP_SOURCE_PARAM = source_cpp="$(EXTNOCYTHONCPP)"
  SETUP = setup_nocython
endif

# C++ compiler flags
EXTCXXFLAGS = $(glucat_extra_cxxflags_pyclical) $(CXXFLAGS)
# set the include path found by configure
EXTAM_CPPFLAGS = $(all_includes)

pyexec_DATA = $(ext_name).*.so

CHECK_OUT = check.out

all: $(pyexec_DATA)

$(pyexec_DATA): $(SETUP_SOURCE_LIST) $(SETUP).py setup_ext.py
	ext_name="$(ext_name)" $(SETUP_SOURCE_PARAM) \
	CXX="$(CXX)" CXXVERSION="$(CXXVERSION)" CFLAGS="" CXXFLAGS="$(EXTCXXFLAGS)" AM_CPPFLAGS="$(EXTAM_CPPFLAGS)" \
	all_includes="$(all_includes)" LDFLAGS="$(USER_LDFLAGS)" LIBRARIES="$(LIBS)" \
	$(PYTHON) $(SETUP).py build_ext --inplace

if use_cython
$(ext_name).cpp: $(pyexec_DATA)

$(EXTNOCYTHONCPP): $(ext_name).cpp
	cp $(ext_name).cpp $@
	$(MAKE) clean

endif

check: $(CHECK_OUT)

$(CHECK_OUT): all
	$(PYTHON) -iq < test.py > $(CHECK_OUT) 2>&1

CLEAN_CYTHON = $(ext_name).cpp $(ext_name).html

CLEANFILES = $(pyexec_DATA) $(CLEAN_CYTHON)

clean-local:
	rm -f *.pyc
	rm -rf build
	rm -rf __pycache__
	rm -f $(CHECK_OUT)

endif
