SHELL := bash

.PHONY: all clean

# Quarto's own project cache (.quarto/) already tracks which entries
# changed since the last render, so this just wraps the render command
# rather than re-implementing Make-level dependency tracking over every
# entries/*.qmd file.
all: _book/index.html

_book/index.html: _quarto.yml index.qmd $(wildcard entries/*.qmd)
	quarto render .

clean:
	rm -rf _book .quarto
