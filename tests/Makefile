
SHELL = /bin/bash -beEu -o pipefail

srcDbLoadGenbankAligns = ../bin/srcDbLoadGenbankAligns
srcDbLoadGenbankSeqs = ../bin/srcDbLoadGenbankSeqs
srcDbLoadGenbankMetadata = ../bin/srcDbLoadGenbankMetadata
srcDbLoadEnsemblAligns = ../bin/srcDbLoadEnsemblAligns
srcDbLoadEnsemblSeqs = ../bin/srcDbLoadEnsemblSeqs
srcDbLoadEnsemblMetadata = ../bin/srcDbLoadEnsemblMetadata
srcDbCheck = ../bin/srcDbCheck
genomeDbLoad = ../bin/genomeDbLoad

genbankTypes = rna refseq est
genbankMetaTypes = rna refseq
ensemblTypes = gencode ensembl

# this obtains the test db to use; gencode must be special cased
#  $(call srcDbLoadTestDb,setName)
srcDbLoadTestDb = $(if $(filter $(1),gencode),hg19,felCat5)


diff = diff

test: srcDbLoadGenbankAlignsTests srcDbLoadGenbankMetadataTests srcDbLoadGenbankSeqsTests \
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
# Run this to update test data from basic tests.  We avoid duplicating tables, as
# metadata and seqs both have align data
##
update: ${genbankTypes:%=%_genbankAlignsUpdate} ${genbankTypes:%=%_genbankSeqsUpdate} ${genbankMetaTypes:%=%_genbankMetadataUpdate} \
	${ensemblTypes:%=%_ensemblAlignsUpdate} ${ensemblTypes:%=%_ensemblSeqsUpdate} ${ensemblTypes:%=%_ensemblMetadataUpdate}

%_genbankAlignsUpdate:
	${MAKE} $*_srcDbLoadGenbankAlignsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadGenbankAlignsTest.db .dump > input/$*_align.sql

%_genbankSeqsUpdate:
	${MAKE} $*_srcDbLoadGenbankSeqsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadGenbankSeqsTest.db '.dump srcSeq' > input/$*_seq.sql

%_genbankMetadataUpdate:
	${MAKE} $*_srcDbLoadGenbankMetadataTest diff=true
	sqlite3 -batch output/$*_srcDbLoadGenbankMetadataTest.db '.dump srcMetadata' > input/$*_metadata.sql

%_ensemblAlignsUpdate:
	${MAKE} $*_srcDbLoadEnsemblAlignsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadEnsemblAlignsTest.db .dump > input/$*_align.sql

%_ensemblSeqsUpdate:
	${MAKE} $*_srcDbLoadEnsemblSeqsTest diff=true
	sqlite3 -batch output/$*_srcDbLoadEnsemblSeqsTest.db '.dump srcSeq' > input/$*_seq.sql

%_ensemblMetadataUpdate:
	${MAKE} $*_srcDbLoadEnsemblMetadataTest diff=true
	sqlite3 -batch output/$*_srcDbLoadEnsemblMetadataTest.db '.dump srcMetadata' > input/$*_metadata.sql

##
# genomeDbLoadTests
##
tree1 = /hive/data/genomes/hg38/bed/multiz20way/hg38.20way.nh
genomeDbLoadTests: genomeDbLoadTest

# felCat3 is projection ensembl
genomeDbLoadTest: mkdirs
	rm -f output/$@.db
	${genomeDbLoad} --hgDb=mm10 --hgDb=hg38 --hgDb=felCat3 --hgDb=canFam3 --tree=${tree1} hgcentraltest output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select * from chains order by srcDb, destDb, chainType' >output/$@.chainsinfo
	sqlite3 -batch output/$@.db '.header on' 'select * from genomeAsms order by hgDb' >output/$@.genomeasms
	${diff} expected/$@.chainsinfo output/$@.chainsinfo
	${diff} expected/$@.genomeasms output/$@.genomeasms

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
