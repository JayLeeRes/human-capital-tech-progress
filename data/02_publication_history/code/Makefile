SHELL := bash

.PHONY: all abstracts clean

# The three core steps below all read AND write ../output/publication_history.csv
# in place (each adds columns to the same file), so plain file-based `make`
# dependencies can't tell them apart by output filename. We use empty
# ".stamp_*" marker files instead, one per step, each depending on the
# previous step's stamp -- this enforces the right order and lets `make`
# skip steps whose stamp is already newer than their own script + inputs.

all: .stamp_own_affiliation

# Pull stage 01's output in as a symlinked input, building stage 01 first
# if needed. Mirrors doctor-tasks' `input/%: ../tasks/X/output/% | input/`
# symlink-recipe convention.
../input/t6_authors.csv: | ../input/
	$(MAKE) -C ../../01_download_journal/code
	ln -sf ../../01_download_journal/output/t6_authors.csv $@

.stamp_fetch: fetch_publication_history.py ../input/t6_authors.csv
	python3 fetch_publication_history.py
	touch $@

.stamp_deep_impact: add_deep_impact_flag.py ../input/deep_impact_journal_list.qmd .stamp_fetch
	python3 add_deep_impact_flag.py
	touch $@

.stamp_own_affiliation: add_own_affiliation.py .stamp_deep_impact
	python3 add_own_affiliation.py
	touch $@

# Optional enrichment step. Depends on stage 03's citation_10y.csv (for its
# eid list), so it must be run AFTER `make -C ../../03_citation_10y/code`,
# not as part of `all` above. Not consumed by anything in stage 04.
abstracts: ../output/abstracts.csv

../output/abstracts.csv: fetch_abstracts.py ../../03_citation_10y/output/citation_10y.csv
	python3 fetch_abstracts.py

clean:
	rm -f .stamp_fetch .stamp_deep_impact .stamp_own_affiliation
	rm -f ../output/publication_history.csv ../output/abstracts.csv
