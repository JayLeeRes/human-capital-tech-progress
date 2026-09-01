SHELL=bash

all: _book/index.html

_book/index.html: _quarto.yml index.qmd $(wildcard entries/*.qmd)
	quarto render .

clean:
	rm -rf _book .quarto
