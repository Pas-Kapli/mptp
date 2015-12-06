# Makefile for creating PDF manual from man file

dist_man_MANS = mptp.1

if HAVE_PS2PDF

doc_DATA = mptp_manual.pdf

mptp_manual.pdf : mptp.1
	TEMP=$$(mktemp temp.XXXXXXXX) ; \
	if [ $$(uname) == "Darwin" ] ; then \
		${SED} -e 's/\\-/-/g' $< | \
		iconv -f UTF-8 -t ISO-8859-1 > $$TEMP ; \
	else \
		${SED} -e 's/\\-/-/g' $< > $$TEMP ; \
	fi ; \
	man -t ./$$TEMP | ps2pdf -sPAPERSIZE=a4 - $@ ; \
	rm $$TEMP

CLEANFILES=mptp_manual.pdf

endif

