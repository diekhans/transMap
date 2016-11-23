
progs = srcDbLoadGenbankAligns srcDbLoadGenbankMetadata srcDbLoadGenbankSeqs \
	srcDbLoadEnsemblAligns srcDbLoadEnsemblSeqs srcDbCheck \
	genomeDbLoad mappingChainBuild mappingChainIndex transMapJob

# getEnsemblMetadata
# FIXME: change this to do all once cleaned out
libs = __init__.py transMapConf.py genomeData.py srcData.py genbank.py ensembl.py \
	phyloTrees.py chainsFinder.py genbankConf.py mappingChainData.py bigTransMap.py

all:

test:
	(cd tests && ${MAKE} test)

lint:
	flake8 ${libs:%=lib/py/transMap/%} ${progs:%=bin/%}

clean:
	(cd tests && ${MAKE} clean)
	rm -f lib/py/transMap/*.pyc
