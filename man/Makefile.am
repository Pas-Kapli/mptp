# Makefile for creating PDF manual from man file

dist_man_MANS = delimit.1

if HAVE_PS2PDF

doc_DATA = delimit_manual.pdf

delimit_manual.pdf : delimit.1
	TEMP=$$(mktemp temp.XXXXXXXX) ; \
	if [ $$(uname) == "Darwin" ] ; then \
		${SED} -e 's/\\-/-/g' $< | \
		iconv -f UTF-8 -t ISO-8859-1 > $$TEMP ; \
	else \
		${SED} -e 's/\\-/-/g' $< > $$TEMP ; \
	fi ; \
	man -t ./$$TEMP | ps2pdf -sPAPERSIZE=a4 - $@ ; \
	rm $$TEMP

CLEANFILES=delimit_manual.pdf

endif

