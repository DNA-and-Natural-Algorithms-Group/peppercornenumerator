tests:
	nosetests -v --with-coverage --cover-package=enumerator --cover-html

profile:
	nosetests -v --with-profile --profile-stats-file stats.pf
