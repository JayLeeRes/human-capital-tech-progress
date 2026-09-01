SHELL=bash

all: paper.pdf

# paper.qmd's numbers are computed by data/04_interactive_fixed_effects's
# pipeline; unlike doctor-tasks' \input{}-based LaTeX files, paper.qmd
# doesn't read those CSVs as literal file inputs at render time (the
# figures are written into the prose by hand), so there's no file-level
# dependency to declare here -- run `make -C ../data/04_interactive_fixed_effects/code`
# by hand first if you want the underlying numbers refreshed.
paper.pdf: paper.qmd ../bib/references.bib ../bib/aer.bst
	quarto render paper.qmd --to pdf

clean:
	rm -f paper.pdf *.aux *.log *.out
