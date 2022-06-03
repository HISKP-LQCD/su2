# Makefile

# file paths and variables
STAGG = ./staggered_fermions.Rmd
INPUT_FILE = ./input_file.Rmd

define generate_html
	Rscript -e 'library(rmarkdown); rmarkdown::render("$(1)", "html_document")'
endef

all: staggered input_file
	make clean


input_file: $(INPUT_FILE)
	$(call generate_html,$(INPUT_FILE))

staggered: $(STAGG)
	$(call generate_html,$(STAGG))

#Rscript -e 'library(rmarkdown); rmarkdown::render("$(STAGG)", "html_document")'


clean:
	@rm *~ *.gz *.fls *.fdb_latexmk *.log *.aux *.xml \
	*.bbl *.blg *.out diagram.* *.bcf *.mp  *.dvi \
	2> /dev/null || true
	

