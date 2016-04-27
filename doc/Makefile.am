AUTOMAKE_OPTIONS = foreign

htmldir = $(docdir)/html
pdfdir = $(docdir)/pdf

if HAVE_DOXYGEN
  
  html_DATA = api/html/*
  pdf_DATA = api/GluCat-API-reference-manual.pdf
  
  glucatHEADERS = ../glucat/*.h
  testHEADERS = ../test/*.h
  PyClicalSOURCE = ../pyclical/glucat.pxd ../pyclical/PyClical.h  ../pyclical/PyClical.pxd ../pyclical/PyClical.pyx
  
doc: html pdf

html: $(html_DATA) 

$(html_DATA): .doxygen.out

.doxygen.out: $(glucatHEADERS) $(testHEADERS) $(PyClicalSOURCE)
	mkdir -p api/html
	mkdir -p api/latex
	$(DOXYGEN) $(DOXYFLAGS) $(doxydir)/Doxyfile > .doxygen.out

pdf: $(pdf_DATA)

$(pdf_DATA): .doxygen.out
	$(MAKE) -C api/latex > api/latex/make.out
	cp api/latex/refman.pdf $@
	rm -rf api/latex

clean-local:
	rm -rf api
	rm -f .doxygen.out

.PHONY: doc html pdf

endif
