
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
transMapStats = ../bin/transMapStats
transMapStatsQCPlot = ../bin/transMapStatsQCPlot
transMapDetails = ../bin/transMapDetails
transMapDetailsCmp = ../bin/transMapDetailsCmp
luigiTransMap = ../bin/luigiTransMap
linkToGbdb = ../bin/linkToGbdb
mkPushList = ../bin/mkPushList

bigTransMapAs = ${HOME}/kent/src/hg/lib/bigTransMap.as

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
	srcDbCheckTests genomeDbLoadTests submitTests


###
# srcDbLoadGenbankAligns
###
srcDbLoadGenbankAlignsTests: ${genbankTypes:%=%_srcDbLoadGenbankAlignsTest}

%_srcDbLoadGenbankAlignsTest: mkdirs
	rm -f output/$@.db
	${srcDbLoadGenbankAligns} --testAccFile=input/$*.acc $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db '.header on' 'select * from srcXRef order by srcAlignId' >output/$@.xref
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
	sqlite3 -batch output/$@.db '.header on' 'select * from srcMetadata order by srcId' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# srcDbLoadEnsemblAligns
##
srcDbLoadEnsemblAlignsTests: ${ensemblTypes:%=%_srcDbLoadEnsemblAlignsTest}

%_srcDbLoadEnsemblAlignsTest: mkdirs
	rm -f output/$@.db
	${srcDbLoadEnsemblAligns} --testAccvFile=input/$*.accv $(call srcDbLoadTestDb,$*) ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db '.header on' 'select * from srcXRef order by srcAlignId' >output/$@.xref
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
	sqlite3 -batch output/$@.db '.header on'  'select * from srcMetadata order by srcId' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# compatiblity with older Ensembls that don't have version numbers; we fake it
srcDbLoadEnsemblCompatTest: mkdirs
	rm -f output/$@.db
	${srcDbLoadEnsemblAligns} --testAccvFile=input/fr3.acc fr3 ensembl output/$@.db
	${srcDbLoadEnsemblSeqs} fr3 ensembl output/$@.db
	${srcDbLoadEnsemblMetadata} --testLimit fr3 ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db '.header on' 'select * from srcXRef order by srcAlignId' >output/$@.xref
	sqlite3 -batch output/$@.db '.header on' 'select * from srcMetadata order by srcId' >output/$@.meta
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
	dataRootDir=output/$@.root ${genomeDbLoad} --hgDb=mm10 --hgDb=hg19 --hgDb=hg38 --hgDb=felCat3 --hgDb=canFam3 input/testTransMapConf.py
	sqlite3 -batch output/$@.root/data/genome.db '.header on' 'select * from chains order by srcHgDb, destHgDb, chainType' >output/$@.chainsinfo
	sqlite3 -batch output/$@.root/data/genome.db '.header on' 'select * from genomeAsms order by hgDb' >output/$@.genomeasms
	${diff} expected/$@.chainsinfo output/$@.chainsinfo
	${diff} expected/$@.genomeasms output/$@.genomeasms

#########################################################################################
mappingTests: mappingChainsBuildSynTest mappingChainsBuildAllTest \
	transMapJobTest transMapBatchTest transMapBuildBigPslTest \
	transMapInfoTest transMapStatsTest transMapDetailsTest \
	transMapStatsQCPlotTest transMapDetailsCmpTest

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
rm -rf output/$@.root
mkdir -p output/$@.root/data/src/hg19 output/$@.root/data/chains/hg19/mm10
cat input/refseq_align.sql input/refseq_seq.sql input/refseq_metadata.sql | sqlite3 -batch output/$@.root/data/src/hg19/hg19.refseq.src.db
sqlite3 -batch output/$@.root/data/genome.db < input/mappingGenomeDb.sql
sqlite3 -batch output/$@.root/data/chains/hg19/mm10/hg19.mm10.all.chain.db < input/mappingChains.sql
cp -f input/mm10.hg19.chain output/$@.root/data/chains/hg19/mm10/hg19.mm10.all.chain
endef

# setup hg19->mm10 test bigPsl
define transMapBigPslSetup
mkdir -p output/$@.root/build/results/mm10/
cp -f expected/transMapBatchTest.preBigPsl output/$@.root/build/results/mm10/mm10.hg19.refseq.preBigPsl
dataRootDir=output/$@.root ${transMapBuildBigPsl} input/testTransMapConf.py mm10 refseq hg19
endef

transMapJobTest: mkdirs
	${transMapEnvSetup}
	dataRootDir=output/$@.root ${transMapJob} input/testTransMapConf.py hg19 mm10 refseq all 1 11 output/transMapJobTest.build/hg19/mm10/refseq/results/parts/1.11.preBigBed
	${diff} expected/$@.preBigPsl output/transMapJobTest.build/hg19/mm10/refseq/results/parts/1.11.preBigBed
	bedToBigBed -type=bed12+20 -as=${bigTransMapAs} -tab output/transMapJobTest.build/hg19/mm10/refseq/results/parts/1.11.preBigBed input/mm10.chrom.sizes output/$@.bigPsl

transMapBatchTest: mkdirs
	${transMapEnvSetup}
	dataRootDir=output/$@.root ${transMapBatch} --directRun input/testTransMapConf.py hg19 mm10 refseq all
	${diff} expected/$@.preBigPsl output/transMapBatchTest.root/build/results/mm10/mm10.hg19.refseq.preBigPsl

transMapBuildBigPslTest: mkdirs
	${transMapBigPslSetup}
	bigBedToBed output/transMapBuildBigPslTest.root/results/mapped/mm10/mm10.refseq.transMapV1.bigPsl output/$@.preBigPsl
	${diff} expected/transMapBatchTest.preBigPsl output/$@.preBigPsl

transMapInfoTest: mkdirs
	rm -rf output/$@.root/data
	mkdir -p output/$@.root/data
	sqlite3 -batch output/$@.root/data/genome.db < input/mappingGenomeDb.sql
	dataRootDir=output/$@.root ${transMapInfo} input/testTransMapConf.py output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

transMapStatsTest: mkdirs
	${transMapEnvSetup}
	sqlite3 -batch output/$@.root/data/genome.db \
	   'update genomeAsms set annotationTypeSet="refseq" where hgDb="hg19";' \
	   'update genomeAsms set annotationTypeSet="" where hgDb="mm10";'
	${transMapBigPslSetup}
	dataRootDir=output/$@.root ${transMapStats} input/testTransMapConf.py output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

transMapDetailsTest: mkdirs
	${transMapEnvSetup}
	${transMapBigPslSetup}
	dataRootDir=output/$@.root ${transMapDetails} input/testTransMapConf.py mm10 output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

transMapStatsQCPlotTest: mkdirs
	${transMapStatsQCPlot} --countChangeTsv=output/$@.change.tsv input/prevOverallStats.tsv input/prevPush.lst input/newOverallStats.tsv input/newPush.lst output/$@.pdf
	${diff} expected/$@.change.tsv output/$@.change.tsv

transMapDetailsCmpTest: mkdirs
	${transMapDetailsCmp} input/transMapDetailsCmp/V5 input/transMapDetailsCmp/V6 output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

#########################################################################################
submitTests: linkToGbdbTest mkPushListTest

submitDestDbs = mm10 hg19
submitBigPslTypes = refseq.transMapV1.bigPsl est.transMapV1.bigPsl \
	rna.transMapV1.bigPsl ensembl.transMapV1.bigPsl \

linkTestBigPsls = $(foreach db,${submitDestDbs},$(foreach tp,${submitBigPslTypes}, mapped/${db}/${db}.${tp}))

linkToGbdbTest: mkdirs
	@rm -rf output/$@.root/data output/$@.root/results output/$@.gbdb
	@mkdir -p output/$@.root/data output/$@.root/results output/$@.gbdb
	for bp in ${linkTestBigPsls:%=output/$@.root/results/%} ; do mkdir -p $$(dirname $${bp}) && touch $${bp} ; done
	sqlite3 -batch output/$@.root/data/genome.db < input/mappingGenomeDb.sql
	dataRootDir=output/$@.root gbdbDir=output/$@.gbdb ${linkToGbdb} input/testTransMapConf.py
	find output/linkToGbdbTest.gbdb -type l | sort > output/$@.files
	${diff} expected/$@.files output/$@.files


mkPushBigPsls = $(foreach db,${submitDestDbs},$(foreach tp,${submitBigPslTypes}, ${db}/transMap/V1/${db}.${tp}))
mkPushListTest: mkdirs
	@rm -rf output/$@.root output/$@.gbdb
	mkdir -p output/$@.root/data output/$@.root/results
	for bp in ${mkPushBigPsls:%=output/$@.gbdb/%} ; do mkdir -p $$(dirname $${bp}) && touch $${bp} ; done
	sqlite3 -batch output/$@.root/data/genome.db < input/mappingGenomeDb.sql
	dataRootDir=output/$@.root gbdbDir=output/$@.gbdb ${mkPushList} --forcePush input/testTransMapConf.py output/$@.pushDbs.lst output/$@.pushGbdbs.lst
	${diff} expected/$@.pushDbs.lst output/$@.pushDbs.lst
	${diff} expected/$@.pushGbdbs.lst output/$@.pushGbdbs.lst

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
	nice ${luigiTransMap} --logConfFile=input/testLuigiLogging.conf --destHgDb=canFam3 --srcHgDb=felCat8 input/testLuigiTransMapConf.py

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
	dataRootDir=output/$@.root ${genomeDbLoad} --hgDb=hg19 --hgDb=mm10 input/testTransMapConf.py
	sqlite3 -batch output/$@.root/data/genome.db \
	   'update chains set chainFile="input/mm10.hg19.chain", netFile="input/mm10.hg19.net" where srcHgDb="hg19"' \
	   '.dump' > input/mappingGenomeDb.sql

# Get overlapping chains in case of new data.  We only use chr7 to keep size down
hg19Sizes = /hive/data/genomes/hg19/chrom.sizes
mm10Sizes = /hive/data/genomes/mm10/chrom.sizes
axtChainHg19Dir = /hive/data/genomes/hg19/bed/lastz.mm10/axtChain
mappingChainsUpdate: ${genbankTypes:%=%_genbankAlignsUpdate} ${gencode_ensemblAlignsUpdate}
	for db in ${genbankTypes:%=output/%_srcDbLoadGenbankAlignsTest.db} output/gencode_srcDbLoadEnsemblAlignsTest.db ; do \
	    sqlite3 -batch $$db 'select * from srcAlign where tName="chr7"' | cut -f 2- ; \
	done > output/$@.hg19.psl
	chainFilter -t=chr7 ${axtChainHg19Dir}/hg19.mm10.all.chain.gz | overlapSelect output/$@.hg19.psl -inFmt=chain /dev/stdin output/hg19.mm10.chain
	chainSwap output/hg19.mm10.chain input/mm10.hg19.chain
	chainNet input/mm10.hg19.chain ${mm10Sizes} ${hg19Sizes} /dev/stdout /dev/null | netSyntenic /dev/stdin input/mm10.hg19.net
	${MAKE} mappingChainsBuildAllTest diff=true
	sqlite3 -batch output/mappingChainsBuildAllTest.all.chain.db '.dump' > input/mappingChains.sql

#########################################################################################
mkdirs:
	@mkdir -p output

clean:
	rm -rf output
