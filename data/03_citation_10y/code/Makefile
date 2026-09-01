SHELL=bash

all: .stamp_fetch

../input/:
	mkdir $@

../input/publication_history.csv: ../../02_publication_history/output/publication_history.csv | ../input/
	ln -sf $< $@

.stamp_fetch: fetch_citation_10y.py ../input/publication_history.csv
	python3 fetch_citation_10y.py
	touch $@

clean:
	rm -f .stamp_fetch ../output/citation_10y.csv ../output/citation_10y_cache.csv
