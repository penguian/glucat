if make_pyclical
  PYCLICAL=pyclical
endif

TIMING = gfft_test products squaring transforms

TOPSUBDIRS = $(PYCLICAL) $(TIMING)

EXAMPLES = test00 test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 test11 test12 test13 test14 test15 test16

TESTOUT = test00/test00.out test01/test01.out test02/test02.out test03/test03.out test04/test04.out test05/test05.out \
          test06/test06.out test07/test07.out test08/test08.out test09/test09.out test10/test10.out test11/test11.out \
          test12/test12.out test13/test13.out test14/test14.out test15/test15.out test16/test16.out

GLUCAT_CONFIG_HEADER = glucat/config.h

EXTRA_DIST = admin doc glucat m4 test test_runtime AUTHORS COPYING ChangeLog configure.ac.in glucat.lsm INSTALL Makefile.am NEWS README TODO

all: docs-am  config.h
	

all: glucat all-header

all-header: $(TOPSUBDIRS) $(CONFIG_HEADER)
	cp $(CONFIG_HEADER) $(GLUCAT_CONFIG_HEADER)

$(TOPSUBDIRS)::
	test -z $(GLUCAT_CONFIG_HEADER) || rm -f $(GLUCAT_CONFIG_HEADER)
	$(MAKE) -C $@

check: check-recursive
	

check: glucat check-out-header

glucat:  $(GLUCAT_CONFIG_HEADER)

$(GLUCAT_CONFIG_HEADER):
	if expr match "$(MAKEFLAGS)" ".* -j"; then \
	  while ! [ -f $(GLUCAT_CONFIG_HEADER) ]; do \
	    sleep 1; \
	  done \
	else \
	  cp $(CONFIG_HEADER) $(GLUCAT_CONFIG_HEADER); \
	fi

check-out-header: $(TESTOUT) $(CONFIG_HEADER)
	cat $(TESTOUT) > test_runtime/test.out
	cp $(CONFIG_HEADER) $(GLUCAT_CONFIG_HEADER)

$(TESTOUT)::
	test -z $(GLUCAT_CONFIG_HEADER) || rm -f $(GLUCAT_CONFIG_HEADER)
	test -z test_runtime/test.out   || rm -f test_runtime/test.out
	$(MAKE) -C $(@D) check

.PHONY: all-header check-out-header

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
# Without this fix aclocal would use tar with parameters 'chof' for "make dist",
# and that would use the old v7 format for tar files.
# That won't work with the long file names in the doc directory.
am__tar = $${TAR-tar} chf - "$$tardir"

dist-hook:
	cd $(top_distdir) && perl admin/am_edit -padmin
	cd $(top_distdir) && $(MAKE) -f admin/Makefile.common subdirs

clean-local:
	test -z $(GLUCAT_CONFIG_HEADER) || rm -f $(GLUCAT_CONFIG_HEADER)
	for d in $(EXAMPLES); do \
	  $(MAKE) -C $$d clean; \
	done

SUBDIRS=$(TOPSUBDIRS)
