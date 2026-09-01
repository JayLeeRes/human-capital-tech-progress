SHELL=bash

all: .stamp_own_affiliation

../input/:
	mkdir $@

../input/t6_authors.csv: ../../01_download_journal/output/t6_authors.csv | ../input/
	ln -sf $< $@

.stamp_fetch: fetch_publication_history.py ../input/t6_authors.csv
	python3 fetch_publication_history.py
	touch $@

.stamp_deep_impact: add_deep_impact_flag.py ../input/deep_impact_journal_list.qmd .stamp_fetch
	python3 add_deep_impact_flag.py
	touch $@

.stamp_own_affiliation: add_own_affiliation.py .stamp_deep_impact
	python3 add_own_affiliation.py
	touch $@

# Optional enrichment step -- not consumed by anything downstream. Needs
# stage 03's citation_10y.csv already built (../$< would go one level
# too far here since input/ and code/ are siblings, not nested).
abstracts: ../output/abstracts.csv

../output/abstracts.csv: fetch_abstracts.py ../../03_citation_10y/output/citation_10y.csv
	python3 fetch_abstracts.py

clean:
	rm -f .stamp_fetch .stamp_deep_impact .stamp_own_affiliation
	rm -f ../output/publication_history.csv ../output/abstracts.csv
