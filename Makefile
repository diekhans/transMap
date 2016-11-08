
progs = getGenbankAligns getGenbankMetadata getGenbankSeqs \
	getEnsemblAligns getEnsemblSeqs  \
	checkTransMapSrcDb

# getEnsemblMetadata
# FIXME: change this to do all once cleaned out
libs = __init__.py genomeDefs.py transMapLite.py genbank.py ensembl.py

all:

test:
	(cd tests && ${MAKE} test)

lint:
	flake8 ${libs:%=lib/py/transMap/%} ${progs:%=bin/%}
