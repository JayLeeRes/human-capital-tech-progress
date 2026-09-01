SHELL := bash

.PHONY: all clean

all: paper.pdf

# paper.qmd's tables/figures are computed from data/04_interactive_fixed_effects's
# output/; make sure that stage is up to date before rendering.
paper.pdf: paper.qmd ../bib/references.bib ../bib/aer.bst
	$(MAKE) -C ../data/04_interactive_fixed_effects/code
	quarto render paper.qmd --to pdf

clean:
	rm -f paper.pdf *.aux *.log *.out
