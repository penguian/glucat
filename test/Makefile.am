####### kdevelop will overwrite this part!!! (begin)##########

EXTRA_DIST = driver.h Makefile.am

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
