# Makefile para pozo_potencial_rectangular.tex
#

FILES = pozo_potencial_rectangular.tex\
	texto/ecuacion_schrodinger.tex\
	texto/estados_ligados.tex\
	texto/estados_no_ligados.tex\
	texto/solucion_numerica.tex

%.pdf:	%.tex $(FILES)
	pdflatex --enable-write18 $<
	pdflatex --enable-write18 $<

all: pozo_potencial_rectangular.pdf

.PHONY: all clean

clean:
	rm -rf *.pdf *.ps *.dvi *.aux *.log *.toc *.gnuplot *.table texto/*.aux

