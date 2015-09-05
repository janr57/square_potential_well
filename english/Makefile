# Makefile for square_potential_well.tex
#

FILES = text/schrodinger_equation.tex

%.pdf:	%.tex $(FILES)
	pdflatex --enable-write18 $<
	pdflatex --enable-write18 $<

all: square_potential_well.pdf

.PHONY: all clean

clean:
	rm -rf *.pdf *.ps *.dvi *.aux *.log *.toc text/*.aux

