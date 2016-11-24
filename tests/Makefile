
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

genbankTypes = rna refseq est
genbankMetaTypes = rna refseq
ensemblTypes = gencode ensembl

# this obtains the test db to use; ensemble must be special cased
#  $(call srcDbLoadTestDb,setName)
srcDbLoadTestDb = $(if $(filter $(1),ensembl),felCat5,hg19)

treeHg38 = /hive/data/genomes/hg38/bed/multiz20way/hg38.20way.nh
treeHg19 = /hive/data/genomes/hg19/bed/multiz4way/44way.nh


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
	${srcDbLoadEnsemblAligns} --testAccvFile=input/$*.accv $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRef' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

##
# srcDbLoadEnsemblSeqs
##
srcDbLoadEnsemblSeqsTests: ${ensemblTypes:%=%_srcDbLoadEnsemblSeqsTest}

%_srcDbLoadEnsemblSeqsTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${srcDbLoadEnsemblSeqs} $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeq order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids


##
# srcDbLoadEnsemblMetadata
##
srcDbLoadEnsemblMetadataTests: ${ensemblTypes:%=%_srcDbLoadEnsemblMetadataTest}

%_srcDbLoadEnsemblMetadataTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${srcDbLoadEnsemblMetadata} --testLimit $(call srcDbLoadTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# srcDbCheck
##
srcDbCheckTests: ${genbankTypes:%=%_srcDbCheckTest}  ${ensemblTypes:%=%_srcDbCheckTest}

%_srcDbCheckTest: mkdirs
	rm -f output/$@.db
	cat input/$*_align.sql input/$*_seq.sql input/$*_metadata.sql | sqlite3 -batch output/$@.db
	${srcDbCheck} $(call srcDbLoadTestDb,$*) $* output/$@.db

##
# genomeDbLoadTests
##
genomeDbLoadTests: genomeDbLoadTest

# felCat3 is projection ensembl
genomeDbLoadTest: mkdirs
	rm -f output/$@.db
	${genomeDbLoad} --hgDb=mm10 --hgDb=hg19 --hgDb=hg38 --hgDb=felCat3 --hgDb=canFam3 --tree=${treeHg38} --tree=${treeHg19} hgcentraltest output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select * from chains order by srcDb, destDb, chainType' >output/$@.chainsinfo
	sqlite3 -batch output/$@.db '.header on' 'select * from genomeAsms order by hgDb' >output/$@.genomeasms
	${diff} expected/$@.chainsinfo output/$@.chainsinfo
	${diff} expected/$@.genomeasms output/$@.genomeasms

#########################################################################################
mappingTests: mappingChainsBuildSynTest mappingChainsBuildAllTest \
	transMapJobTest

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

transMapJobTest: mkdirs
	rm -f output/$@.src.db output/$@.chain.db
	cat input/rna_align.sql input/rna_seq.sql input/rna_metadata.sql | sqlite3 -batch output/$@.src.db
	sqlite3 -batch output/$@.chain.db < input/mappingChains.sql
	cp -f input/mm10.hg19.chain output/$@.chain
	${transMapJob} hg19 output/$@.src.db 1 11 output/$@.chain mm10 output/$@.preBigPsl
	${diff} expected/$@.preBigPsl output/$@.preBigPsl
	bedToBigBed -type=bed12+20 -as=../etc/bigTransMap.as -tab output/$@.preBigPsl input/mm10.chrom.sizes output/$@.bigPsl

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
	${genomeDbLoad} --hgDb=hg19 --hgDb=mm10 --tree=${treeHg19} hgcentraltest output/$@.db
	sqlite3 -batch output/$@.db \
	   'update chains set chainFile="input/mm10.hg19.chain", netFile="input/mm10.hg19.net" where srcDb="hg19"' \
	   '.dump' > input/mappingGenomeDb.sql

mappingChainsUpdate: mkdirs
	${MAKE} mappingChainsBuildAllTest diff=true
	sqlite3 -batch output/mappingChainsBuildAllTest.all.chain.db '.dump' > input/mappingChains.sql

#########################################################################################
mkdirs:
	@mkdir -p output

clean:
	rm -rf output
