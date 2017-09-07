.PHONY: tests
.PHONY: profile
.PHONY: docs
.PHONY: README

tests:
	python -m nose -v --with-coverage --cover-package=enumerator --cover-html tests/test_reactions_SB.py tests/test_utils.py tests/test_objects.py

profile:
	nosetests -v --with-profile --profile-stats-file stats.pf

docs: 
	cd docs && $(MAKE) html

README: README.pdf README.html
	

README.pdf: README.md
	pandoc README.md -o README.pdf

README.html: README.md
	pandoc -s README.md -o README.html
