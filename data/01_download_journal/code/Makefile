SHELL=bash

all: .stamp_fetch

.stamp_fetch: fetch_t6_papers.py
	python3 fetch_t6_papers.py
	touch $@

clean:
	rm -f .stamp_fetch ../output/t6_papers.csv ../output/t6_authors.csv
