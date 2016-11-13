
progs = getGenbankAligns getGenbankMetadata getGenbankSeqs \
	getEnsemblAligns getEnsemblSeqs checkTransMapSrcDb \
	getGenomeInfo

# getEnsemblMetadata
# FIXME: change this to do all once cleaned out
libs = __init__.py genomeDefs.py srcData.py genbank.py ensembl.py \
	phyloTrees.py chainsFinder.py

all:

test:
	(cd tests && ${MAKE} test)

lint:
	flake8 ${libs:%=lib/py/transMap/%} ${progs:%=bin/%}
