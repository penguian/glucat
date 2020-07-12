AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = Makefile.in

htmldir = $(docdir)/html
pdfdir = $(docdir)/pdf

if HAVE_DOXYGEN
  
  nodist_html_DATA = api/html/*
  nodist_pdf_DATA = api/GluCat-API-reference-manual.pdf
  
  glucatHEADERS = ../glucat/*.h
  testHEADERS = ../test/*.h
  PyClicalSOURCE = ../pyclical/glucat.pxd ../pyclical/PyClical.h  ../pyclical/PyClical.pxd ../pyclical/PyClical.pyx
  
doc: html pdf

html: $(nodist_html_DATA) 

$(nodist_html_DATA): .doxygen.out

.doxygen.out: $(glucatHEADERS) $(testHEADERS) $(PyClicalSOURCE)
	mkdir -p api/html
	mkdir -p api/latex
	$(DOXYGEN) -u $(doxydir)/Doxyfile > .doxygen-u.out 2>&1
	$(DOXYGEN) $(DOXYFLAGS) $(doxydir)/Doxyfile > .doxygen.out

pdf: $(nodist_pdf_DATA)

$(nodist_pdf_DATA): .doxygen.out
	sed -i "s/operator\%/operator\\\%/g" api/latex/*.tex
	$(MAKE) -C api/latex > api/latex/make.out
	cp api/latex/refman.pdf $@
	rm -rf api/latex

clean-local:
	rm -rf api
	rm -f .doxygen-u.out
	rm -f .doxygen.out

.PHONY: doc html pdf

endif
