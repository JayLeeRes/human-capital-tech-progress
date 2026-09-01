SHELL := bash

.PHONY: all clean

all: .stamp_fetch

# Pull stage 02's output in as a symlinked input, building stage 02 first
# if needed.
../input/publication_history.csv: | ../input/
	$(MAKE) -C ../../02_publication_history/code
	ln -sf ../../02_publication_history/output/publication_history.csv $@

# One script call writes both outputs (the resumable raw cache and the
# final aggregated file). Routed through a stamp file for the same reason
# as stage 01 (GNU Make 3.81 has no grouped multi-target rules).
.stamp_fetch: fetch_citation_10y.py ../input/publication_history.csv
	python3 fetch_citation_10y.py
	touch $@

clean:
	rm -f .stamp_fetch ../output/citation_10y.csv ../output/citation_10y_cache.csv
