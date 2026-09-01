SHELL := bash

.PHONY: all clean

# fetch_t6_papers.py writes both outputs from one Scopus pull. GNU Make 3.81
# (macOS default) doesn't support grouped multi-target rules ("&:"), and a
# plain multi-target rule would run the recipe once per target -- so we
# route both through a single stamp file instead, same as stage 02.

all: .stamp_fetch

.stamp_fetch: fetch_t6_papers.py
	python3 fetch_t6_papers.py
	touch $@

clean:
	rm -f .stamp_fetch ../output/t6_papers.csv ../output/t6_authors.csv
