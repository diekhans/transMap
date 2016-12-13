
SHELL = /bin/bash -beEu -o pipefail

srcDbLoadGenbankAligns = ../bin/srcDbLoadGenbankAligns
srcDbLoadGenbankSeqs = ../bin/srcDbLoadGenbankSeqs
srcDbLoadGenbankMetadata = ../bin/srcDbLoadGenbankMetadata
srcDbLoadEnsemblAligns = ../bin/srcDbLoadEnsemblAligns
srcDbLoadEnsemblSeqs = ../bin/srcDbLoadEnsemblSeqs
srcDbLoadEnsemblMetadata = ../bin/srcDbLoadEnsemblMetadata
srcDbCheck = ../bin/srcDbCheck
genomeDbLoad = ../bin/genomeDbLoad
mappingChainBuild = ../bin/mappingChainBuild
transMapJob = ../bin/transMapJob
transMapBatch = ../bin/transMapBatch
transMapBuildBigPsl = ../bin/transMapBuildBigPsl
transMapInfo = ../bin/transMapInfo
luigiTransMap = ../bin/luigiTransMap

genbankTypes = rna refseq est
genbankMetaTypes = rna refseq
# gencode is a pseudo-type for constructing test cases
ensemblTypes = gencode ensembl

# this obtains the annotaton type, changing gencode to ensembl, which is used
# by test, to ensembl
#  $(call getAnnType,setName)
getAnnType = $(if $(filter $(1),gencode),ensembl,$(1))


# this obtains the test db to use; ensembl must be special cased
#  $(call srcDbLoadTestDb,setName)
srcDbLoadTestDb = $(if $(filter $(1),ensembl),felCat5,hg19)

diff = diff

test: srcDbTests mappingTests

#########################################################################################
srcDbTests: srcDbLoadGenbankAlignsTests srcDbLoadGenbankMetadataTests srcDbLoadGenbankSeqsTests \
	srcDbLoadEnsemblAlignsTests srcDbLoadEnsemblSeqsTests srcDbLoadEnsemblMetadataTests \
	srcDbCheckTests genomeDbLoadTests


###
# srcDbLoadGenbankAligns
###
srcDbLoadGenbankAlignsTests: ${genbankTypes:%=%_srcDbLoadGenbankAlignsTest}

%_srcDbLoadGenbankAlignsTest: mkdirs
	rm -f output/$@.db
	${srcDbLoadGenbankAligns} --testAccFile=input/$*.acc $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRef' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

###
# srcDbLoadGenbankSeqs
###
srcDbLoadGenbankSeqsTests: ${genbankTypes:%=%_srcDbLoadGenbankSeqsTest}

%_srcDbLoadGenbankSeqsTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${srcDbLoadGenbankSeqs} $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeq order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids

###
# srcDbLoadGenbankMetadata
###
srcDbLoadGenbankMetadataTests: ${genbankMetaTypes:%=%_srcDbLoadGenbankMetadataTest}

%_srcDbLoadGenbankMetadataTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${srcDbLoadGenbankMetadata} --testLimit $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# srcDbLoadEnsemblAligns
##
srcDbLoadEnsemblAlignsTests: ${ensemblTypes:%=%_srcDbLoadEnsemblAlignsTest}

%_srcDbLoadEnsemblAlignsTest: mkdirs
	rm -f output/$@.db
	${srcDbLoadEnsemblAligns} --testAccvFile=input/$*.accv $(call srcDbLoadTestDb,$*) ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRef' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

##
# srcDbLoadEnsemblSeqs
##
srcDbLoadEnsemblSeqsTests: ${ensemblTypes:%=%_srcDbLoadEnsemblSeqsTest} \
	srcDbLoadEnsemblCompatTest

%_srcDbLoadEnsemblSeqsTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${srcDbLoadEnsemblSeqs} $(call srcDbLoadTestDb,$*) ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeq order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids


##
# srcDbLoadEnsemblMetadata
##
srcDbLoadEnsemblMetadataTests: ${ensemblTypes:%=%_srcDbLoadEnsemblMetadataTest}

%_srcDbLoadEnsemblMetadataTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${srcDbLoadEnsemblMetadata} --testLimit $(call srcDbLoadTestDb,$*) ensembl output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# compatiblity with older Ensembls that don't have version numbers; we fake it
srcDbLoadEnsemblCompatTest: mkdirs
	rm -f output/$@.db
	${srcDbLoadEnsemblAligns} --testAccvFile=input/fr3.acc fr3 ensembl output/$@.db
	${srcDbLoadEnsemblSeqs} fr3 ensembl output/$@.db
	${srcDbLoadEnsemblMetadata} --testLimit fr3 ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRef' >output/$@.xref
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeq order by name' >output/$@.fa.ids
	${srcDbCheck} fr3 ensembl output/$@.db
	${diff} expected/$@.fa.ids output/$@.fa.ids
	${diff} expected/$@.meta output/$@.meta
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

##
# srcDbCheck
##
srcDbCheckTests: ${genbankTypes:%=%_srcDbCheckTest}  ${ensemblTypes:%=%_srcDbCheckTest}

%_srcDbCheckTest: mkdirs
	rm -f output/$@.db
	cat input/$*_align.sql input/$*_seq.sql input/$*_metadata.sql | sqlite3 -batch output/$@.db
	${srcDbCheck} $(call srcDbLoadTestDb,$*) $(call getAnnType,$*) output/$@.db

##
# genomeDbLoadTests
##
genomeDbLoadTests: genomeDbLoadTest

# felCat3 is projection ensembl
genomeDbLoadTest: mkdirs
	rm -rf output/$@.data
	dataDir=output/$@.data ${genomeDbLoad} --hgDb=mm10 --hgDb=hg19 --hgDb=hg38 --hgDb=felCat3 --hgDb=canFam3 input/testTransMapConf.py
	sqlite3 -batch output/$@.data/genome.db '.header on' 'select * from chains order by srcHgDb, destHgDb, chainType' >output/$@.chainsinfo
	sqlite3 -batch output/$@.data/genome.db '.header on' 'select * from genomeAsms order by hgDb' >output/$@.genomeasms
	${diff} expected/$@.chainsinfo output/$@.chainsinfo
	${diff} expected/$@.genomeasms output/$@.genomeasms

#########################################################################################
mappingTests: mappingChainsBuildSynTest mappingChainsBuildAllTest \
	transMapJobTest transMapBatchTest transMapBuildBigPslTest \
	transMapInfoTest

mappingChainsBuildSynTest: mkdirs
	rm -f output/$@.genome.db
	sqlite3 -batch output/$@.genome.db < input/mappingGenomeDb.sql
	${mappingChainBuild} output/$@.genome.db hg19 mm10 syn output/$@.syn.chain
	sqlite3 -batch output/$@.syn.chain.db  'select * from mappingChains' >output/$@.chainidx
	${diff} expected/$@.chainidx output/$@.chainidx

# never use `all' in actually transmap, but excerise rbest code path
mappingChainsBuildAllTest: mkdirs
	rm -f output/$@.genome.db
	sqlite3 -batch output/$@.genome.db < input/mappingGenomeDb.sql
	${mappingChainBuild} output/$@.genome.db hg19 mm10 all output/$@.all.chain
	sqlite3 -batch output/$@.all.chain.db  'select * from mappingChains' >output/$@.chainidx
	${diff} expected/$@.chainidx output/$@.chainidx

# sets up test files for batch/job tests
define transMapEnvSetup
rm -rf output/$@.data output/$@.build
mkdir -p output/$@.data/src/hg19 output/$@.build/chains/mm10/hg19
cat input/rna_align.sql input/rna_seq.sql input/rna_metadata.sql | sqlite3 -batch output/$@.data/src/hg19/hg19.rna.src.db
sqlite3 -batch output/$@.build/chains/mm10/hg19/hg19.mm10.all.chain.db < input/mappingChains.sql
cp -f input/mm10.hg19.chain output/$@.build/chains/mm10/hg19/hg19.mm10.all.chain
endef

transMapJobTest: mkdirs
	${transMapEnvSetup}
	buildTmpDir=output/$@.build dataDir=output/$@.data ${transMapJob} input/testTransMapConf.py hg19 mm10 rna all 1 11 output/transMapJobTest.build/hg19/mm10/rna/results/parts/1.11.preBigBed
	${diff} expected/$@.preBigPsl output/transMapJobTest.build/hg19/mm10/rna/results/parts/1.11.preBigBed
	bedToBigBed -type=bed12+20 -as=../etc/bigTransMap.as -tab output/transMapJobTest.build/hg19/mm10/rna/results/parts/1.11.preBigBed input/mm10.chrom.sizes output/$@.bigPsl

transMapBatchTest: mkdirs
	${transMapEnvSetup}
	buildTmpDir=output/$@.build dataDir=output/$@.data ${transMapBatch} --directRun input/testTransMapConf.py hg19 mm10 rna all
	${diff} expected/$@.preBigPsl output/transMapBatchTest.build/results/mm10/mm10.hg19.rna.preBigPsl

transMapBuildBigPslTest: mkdirs
	mkdir -p output/$@.build/results/mm10/
	cp -f expected/transMapBatchTest.preBigPsl output/$@.build/results/mm10/mm10.hg19.rna.preBigPsl
	buildTmpDir=output/$@.build dataDir=output/$@.data ${transMapBuildBigPsl} input/testTransMapConf.py mm10 rna hg19
	bigBedToBed output/transMapBuildBigPslTest.data/mapped/mm10/mm10.rna.transMap.bigPsl output/$@.preBigPsl
	${diff} expected/transMapBatchTest.preBigPsl output/$@.preBigPsl

transMapInfoTest: mkdirs
	rm -rf output/$@.data
	mkdir -p output/$@.data
	sqlite3 -batch output/$@.data/genome.db < input/mappingGenomeDb.sql
	dataDir=output/$@.data ${transMapInfo} input/testTransMapConf.py output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

#########################################################################################

.PHONY: bigTest
bigTest: bigTestGenomeDb bigTestInfo

bigTestGenomeDb: bigTest/data/genome.db

bigTest/data/genome.db:
	dataDir=bigTest ${genomeDbLoad} input/testTransMapConf.py

bigTestInfo: bigTest/data/genome.db
	dataDir=bigTest/data ${transMapInfo} input/testTransMapConf.py bigTest/transmap-info.tsv

bigClean:
	rm -rf bigTest

#########################################################################################

.PHONY: luigiTest
luigiTest:
	nice ${luigiTransMap} --logConfFile=../etc/logging.conf --destHgDb=canFam3 --srcHgDb=felCat8 input/testLuigiTransMapConf.py

luigiBatchClean:
	rm -rf luigiTest/build/*/*/{est,rna,refseq,ensembl,gencode}

luigiClean:
	rm -rf luigiTest

#########################################################################################

##
# Run this to update test data from basic tests.  We avoid duplicating tables, as
# metadata and seqs test both have align data
##
update: ${genbankTypes:%=%_genbankAlignsUpdate} ${genbankTypes:%=%_genbankSeqsUpdate} ${genbankMetaTypes:%=%_genbankMetadataUpdate} \
	${ensemblTypes:%=%_ensemblAlignsUpdate} ${ensemblTypes:%=%_ensemblSeqsUpdate} ${ensemblTypes:%=%_ensemblMetadataUpdate} \
	mappingGenomeDbUpdate mappingChainsUpdate

%_genbankAlignsUpdate:
	${MAKE} $*_srcDbLoadGenbankAlignsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadGenbankAlignsTest.db .dump > input/$*_align.sql

%_genbankSeqsUpdate: %_genbankAlignsUpdate
	${MAKE} $*_srcDbLoadGenbankSeqsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadGenbankSeqsTest.db '.dump srcSeq' > input/$*_seq.sql

%_genbankMetadataUpdate: %_genbankAlignsUpdate
	${MAKE} $*_srcDbLoadGenbankMetadataTest diff=true
	sqlite3 -batch output/$*_srcDbLoadGenbankMetadataTest.db '.dump srcMetadata' > input/$*_metadata.sql

%_ensemblAlignsUpdate:
	${MAKE} $*_srcDbLoadEnsemblAlignsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadEnsemblAlignsTest.db .dump > input/$*_align.sql

%_ensemblSeqsUpdate: %_ensemblAlignsUpdate
	${MAKE} $*_srcDbLoadEnsemblSeqsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadEnsemblSeqsTest.db '.dump srcSeq' > input/$*_seq.sql

%_ensemblMetadataUpdate: %_ensemblAlignsUpdate
	${MAKE} $*_srcDbLoadEnsemblMetadataTest diff=true
	sqlite3 -batch output/$*_srcDbLoadEnsemblMetadataTest.db '.dump srcMetadata' > input/$*_metadata.sql

mappingGenomeDbUpdate: mkdirs
	dataDir=output/$@.data ${genomeDbLoad} --hgDb=hg19 --hgDb=mm10 input/testTransMapConf.py
	sqlite3 -batch output/$@.data/genome.db \
	   'update chains set chainFile="input/mm10.hg19.chain", netFile="input/mm10.hg19.net" where srcHgDb="hg19"' \
	   '.dump' > input/mappingGenomeDb.sql

mappingChainsUpdate: mkdirs
	${MAKE} mappingChainsBuildAllTest diff=true
	sqlite3 -batch output/mappingChainsBuildAllTest.all.chain.db '.dump' > input/mappingChains.sql

#########################################################################################
mkdirs:
	@mkdir -p output

clean:
	rm -rf output
