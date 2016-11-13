
SHELL = /bin/bash -beEu -o pipefail

getGenbankAligns = ../../bin/getGenbankAligns
getGenbankSeqs = ../../bin/getGenbankSeqs
getGenbankMetadata = ../../bin/getGenbankMetadata
getEnsemblAligns = ../../bin/getEnsemblAligns
getEnsemblSeqs = ../../bin/getEnsemblSeqs
getEnsemblMetadata = ../../bin/getEnsemblMetadata
checkTransMapSrcDb = ../../bin/checkTransMapSrcDb
getGenomeInfo = ../../bin/getGenomeInfo

genbankTypes = rna refseq est
genbankMetaTypes = rna refseq
ensemblTypes = gencode ensembl

# this obtains the test db to use; gencode must be special cased
#  $(call getTestDb,setName)
getTestDb = $(if $(filter $(1),gencode),hg19,felCat5)


diff = diff

test: getGenbankAlignsTests getGenbankMetadataTests getGenbankSeqsTests \
	getEnsemblAlignsTests getEnsemblSeqsTests getEnsemblMetadataTests \
	checkTransMapSrcDbTests getGenomeInfoTests


###
# getGenbankAligns
###
getGenbankAlignsTests: ${genbankTypes:%=%_getGenbankAlignsTest}

%_getGenbankAlignsTest: mkdirs
	rm -f output/$@.db
	${getGenbankAligns} --limit=100 $(call getTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRef' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

###
# getGenbankSeqs
###
getGenbankSeqsTests: ${genbankTypes:%=%_getGenbankSeqsTest}

%_getGenbankSeqsTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${getGenbankSeqs} $(call getTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeq order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids

###
# getGenbankMetadata
###
getGenbankMetadataTests: ${genbankMetaTypes:%=%_getGenbankMetadataTest}

%_getGenbankMetadataTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${getGenbankMetadata} --testLimit $(call getTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# getEnsemblAligns
##
getEnsemblAlignsTests: ${ensemblTypes:%=%_getEnsemblAlignsTest}

%_getEnsemblAlignsTest: mkdirs
	rm -f output/$@.db
	${getEnsemblAligns} --limit=100 $(call getTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAlign' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRef' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

##
# getEnsemblSeqs
##
getEnsemblSeqsTests: ${ensemblTypes:%=%_getEnsemblSeqsTest}

%_getEnsemblSeqsTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${getEnsemblSeqs} $(call getTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeq order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids


##
# getEnsemblMetadata
##
getEnsemblMetadataTests: ${ensemblTypes:%=%_getEnsemblMetadataTest}

%_getEnsemblMetadataTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/$*_align.sql
	${getEnsemblMetadata} --testLimit $(call getTestDb,$*) $* output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# checkTransMapSrcDb
##
checkTransMapSrcDbTests: ${genbankTypes:%=%_checkTransMapSrcDbTest}  ${ensemblTypes:%=%_checkTransMapSrcDbTest}

%_checkTransMapSrcDbTest: mkdirs
	rm -f output/$@.db
	cat input/$*_align.sql input/$*_seq.sql input/$*_metadata.sql | sqlite3 -batch output/$@.db
	${checkTransMapSrcDb} $(call getTestDb,$*) $* output/$@.db


##
# Run this to update test data from basic tests.  We avoid duplicating tables, as
# metadata and seqs both have align data
##
update: ${genbankTypes:%=%_genbankAlignsUpdate} ${genbankTypes:%=%_genbankSeqsUpdate} ${genbankMetaTypes:%=%_genbankMetadataUpdate} \
	${ensemblTypes:%=%_ensemblAlignsUpdate} ${ensemblTypes:%=%_ensemblSeqsUpdate} ${ensemblTypes:%=%_ensemblMetadataUpdate}

%_genbankAlignsUpdate:
	${MAKE} $*_getGenbankAlignsTest diff=true
	sqlite3 -batch output/$*_getGenbankAlignsTest.db .dump > input/$*_align.sql

%_genbankSeqsUpdate:
	${MAKE} $*_getGenbankSeqsTest diff=true
	sqlite3 -batch output/$*_getGenbankSeqsTest.db '.dump srcSeq' > input/$*_seq.sql

%_genbankMetadataUpdate:
	${MAKE} $*_getGenbankMetadataTest diff=true
	sqlite3 -batch output/$*_getGenbankMetadataTest.db '.dump srcMetadata' > input/$*_metadata.sql

%_ensemblAlignsUpdate:
	${MAKE} $*_getEnsemblAlignsTest diff=true
	sqlite3 -batch output/$*_getEnsemblAlignsTest.db .dump > input/$*_align.sql

%_ensemblSeqsUpdate:
	${MAKE} $*_getEnsemblSeqsTest diff=true
	sqlite3 -batch output/$*_getEnsemblSeqsTest.db '.dump srcSeq' > input/$*_seq.sql

%_ensemblMetadataUpdate:
	${MAKE} $*_getEnsemblMetadataTest diff=true
	sqlite3 -batch output/$*_getEnsemblMetadataTest.db '.dump srcMetadata' > input/$*_metadata.sql

##
# getGenomeInfoTests
##
tree1 = /hive/data/genomes/hg38/bed/multiz20way/hg38.20way.nh
getGenomeInfoTests: getGenomeInfoTest

getGenomeInfoTest: mkdirs
	rm -f output/$@.db
	${getGenomeInfo} --hgDb=mm10 --hgDb=hg38 --tree=${tree1} hgcentraltest output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select * from chains order by srcDb, destDb, chainType' >output/$@.chainsinfo
	${diff} expected/$@.chainsinfo output/$@.chainsinfo

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
