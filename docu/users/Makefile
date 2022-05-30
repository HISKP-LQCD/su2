# Makefile

all: staggered
	make clean

STAGG = ./staggered_fermions.Rmd
staggered: $(STAGG)
	Rscript -e 'library(rmarkdown); rmarkdown::render("$(STAGG)", "html_document")'

clean:
	@rm *~ *.gz *.fls *.fdb_latexmk *.log *.aux *.xml \
	*.bbl *.blg *.out diagram.* *.bcf *.mp  *.dvi \
	2> /dev/null || true
	

