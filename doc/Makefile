LATEX = latex
DVIPDF = dvipdf
BIBTEX = bibtex
PDFLATEX = pdflatex

all: delimit

delimit: delimit.tex
	$(LATEX) $+
	$(BIBTEX) $@
	$(DVIPDF) $@
