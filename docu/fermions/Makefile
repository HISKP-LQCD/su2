# Makefile

CC = pdflatex
TEX = main.tex
PDF = main.pdf
BCF = main.bcf
BIB = biblio.bib

all: $(BCF)
	make clean

show: $(PDF) $(TEX)
	evince $(PDF) &

$(BCF): $(BIB) $(TEX)
	make -B $(PDF)
	biber $(BCF)
	make -B $(PDF)
	make -B $(PDF)	

$(PDF): $(TEX)
	$(CC) $(TEX)


clean:
	@rm *~ *.gz *.fls *.fdb_latexmk *.log *.aux *.xml \
	*.bbl *.blg *.out diagram.* *.bcf *.mp  \
	2> /dev/null || true
	

