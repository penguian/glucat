AUTOMAKE_OPTIONS = foreign

htmldir = $(docdir)/html
pdfdir = $(docdir)/pdf

if HAVE_DOXYGEN
  
  html_DATA = api/html/*
  pdf_DATA = api/GluCat-API-reference-manual.pdf
  
html: $(html_DATA)

$(html_DATA): doxygen.out

doxygen.out: $(doxydir)/Doxyfile
	mkdir -p api/html
	mkdir -p api/latex
	$(DOXYGEN) $(DOXYFLAGS) $< > doxygen.out

pdf: $(pdf_DATA)

$(pdf_DATA): html
	$(MAKE) -C api/latex > api/latex/make.out
	cp api/latex/refman.pdf $@
	$(MAKE) clean

  CLEANFILES = doxygen.out

clean-local:
	rm -rf api/latex

distclean-local:
	rm -rf api

.PHONY: html pdf

endif
