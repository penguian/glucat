if make_pyclical
  PYCLICAL=pyclical
endif

TOPSUBDIRS = glucat test transforms gfft_test products squaring test00 test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 test11 test12 test13 test14 test15 test16 $(PYCLICAL)

EXTRA_DIST = admin doc m4 test_runtime AUTHORS COPYING ChangeLog configure.ac.in glucat.lsm INSTALL Makefile.am NEWS README TODO

# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS = foreign

$(top_srcdir)/configure.ac: configure.ac.in $(top_srcdir)/subdirs
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common configure.ac ;

$(top_srcdir)/configure.files:
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common configure.files ;

$(top_srcdir)/subdirs:
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common subdirs ;

MAINTAINERCLEANFILES = subdirs configure.in acinclude.m4 configure.files

cvs:
	$(MAKE) -f admin/Makefile.common cvs

# Fix to "make dist" which needs tar.
# Without this fix aclocal would use tar with parameters 'chof' for "mke dist",
# and that would use the old v7 format for tar files.
# That won't work with the long file names in the doc directory.
am__tar = $${TAR-tar} chf - "$$tardir"

dist-hook:
	cd $(top_distdir) && perl admin/am_edit -padmin
	cd $(top_distdir) && $(MAKE) -f admin/Makefile.common subdirs


SUBDIRS=$(TOPSUBDIRS)
