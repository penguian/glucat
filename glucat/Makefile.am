####### kdevelop will overwrite this part!!! (begin)##########


EXTRA_DIST = clifford_algebra.h clifford_algebra_imp.h errors.h errors_imp.h framed_multi.h framed_multi_imp.h generation.h generation_imp.h global.h glucat.h glucat_imp.h index_set.h index_set_imp.h long_double.h matrix.h matrix_imp.h matrix_multi.h matrix_multi_imp.h portability.h scalar.h 

####### kdevelop will overwrite this part!!! (end)############
AUTOMAKE_OPTIONS = foreign
# set the include path found by configure
includedir	 = ${prefix}/include/glucat
INCLUDES = $(all_includes)
SUBDIRS =
make-install-dirs:
	-if test '!' -d $(includedir); then mkdir -p $(includedir); fi
install-data-hook: make-install-dirs
	-@ echo Installing $(includedir) ; \
	$(INSTALL_DATA) *.h $(includedir)
