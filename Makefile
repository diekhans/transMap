
progs = getGenbankAligns getGenbankMetadata getGenbankSeqs \
	checkTransMapSrcDb
libs = __init__.py genomeDefs.py transMapLite.py

all:

test:
	(cd tests && ${MAKE} test)

lint:
	flake8 ${libs:%=lib/py/transMap/%} ${progs:%=bin/%}
