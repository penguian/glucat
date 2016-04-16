AUTOMAKE_OPTIONS = foreign

if HAVE_DOXYGEN

  PDFDOC = api/GluCat-API-reference-manual.pdf
  doc_pdf_DATA = $(PDFDOC)
  doc_pdfdir = $(docdir)/pdf
  
  HTMLDOC = api/html
  doc_html_DATA = $(HTMLDOC)/*
  doc_htmldir = $(docdir)/html

$(HTMLDOC)/*: doxygen.out

doxygen.out: @doxydir@/Doxyfile
	mkdir -p api/html
	mkdir -p api/latex
	cd @doxydir@ && $(DOXYGEN) $(DOXYFLAGS) $< > doc/doxygen.out

$(PDFDOC): doxygen.out
	make -C api/latex > make-pdf.out
	cp api/latex/refman.pdf $(PDFDOC)

  CLEANFILES = doxygen.out make-pdf.out

clean-local:
	rm -rf api

endif
