# Makefile

# file paths and variables
STAGG = ./staggered_fermions.Rmd
u1 = ./u1.Rmd

define generate_html
	Rscript -e 'library(rmarkdown); rmarkdown::render("$(1)", "html_document")'
endef

all: staggered u1
	make clean

u1: $(u1)
	$(call generate_html,$(u1))

staggered: $(STAGG)
	$(call generate_html,$(STAGG))


clean:
	@rm *~ *.gz *.fls *.fdb_latexmk *.log *.aux *.xml \
	*.bbl *.blg *.out diagram.* *.bcf *.mp  *.dvi \
	2> /dev/null || true
	

