EXTRA_DIST = pyclical_demo.py pyclical_demo.out Makefile.am

AUTOMAKE_OPTIONS = foreign

ext_name = PyClical

# C++ compiler flags
glucat_extra_cxxflags_pyclical = @glucat_extra_cxxflags_pyclical@
EXTCXXFLAGS = $(glucat_extra_cxxflags_pyclical) $(CXXFLAGS)
# set the include path found by configure
INCLUDES = $(all_includes)
# the library search path.
LDFLAGS = $(all_libraries)

pyexec_DATA = $(ext_name).so

all: $(pyexec_DATA) $(ext_name).html

$(ext_name).so: $(ext_name).cpp $(ext_name).h ../glucat/glucat.h ../glucat/glucat_imp.h
	ext_name="$(ext_name)" source_cpp="$(ext_name).cpp" CXXFLAGS="$(EXTCXXFLAGS)" INCLUDES="$(INCLUDES)" LDFLAGS="$(LDFLAGS)" python setup_nocython.py build_ext --inplace

$(ext_name).cpp $(ext_name).html: setup.py $(ext_name).pyx $(ext_name).pxd glucat.pxd $(ext_name).h ../glucat/glucat.h ../glucat/glucat_imp.h
	ext_name="$(ext_name)" source_pyx="$(ext_name).pyx" CXXFLAGS="$(EXTCXXFLAGS)" INCLUDES="$(INCLUDES)" LDFLAGS="$(LDFLAGS)" python setup.py build_ext --inplace

clean:
	rm -f $(ext_name).so $(ext_name).cpp $(ext_name).html
	rm -rf build