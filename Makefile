.PHONY: tests
.PHONY: profile
.PHONY: docs

tests:
	nosetests -v --with-coverage --cover-package=enumerator --cover-html test_condense.py test_enumerator.py test_input.py test_output.py test_reactions.py test_utils.py
#	nosetests -v --with-coverage --cover-package=enumerator --cover-html

profile:
	nosetests -v --with-profile --profile-stats-file stats.pf

docs: 
	cd docs && $(MAKE) html

README.pdf: README.md
	pandoc README.md -o README.pdf

README.html: README.md
	pandoc -s README.md -o README.html
