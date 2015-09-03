EXTRA_DIST = *.h

AUTOMAKE_OPTIONS = foreign
# set the include path found by configure
includedir	 = ${prefix}/include/glucat
AM_CPPFLAGS = $(all_includes)
SUBDIRS =
make-install-dirs:
	-if test '!' -d $(includedir); then mkdir -p $(includedir); fi
install-data-hook: make-install-dirs
	-@ echo Installing $(includedir) ; \
	$(INSTALL_DATA) *.h $(includedir)
