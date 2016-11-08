
getGenbankAligns = ../../bin/getGenbankAligns
getGenbankSeqs = ../../bin/getGenbankSeqs
getGenbankMetadata = ../../bin/getGenbankMetadata
getEnsemblAligns = ../../bin/getEnsemblAligns
getEnsemblSeqs = ../../bin/getEnsemblSeqs
getEnsemblMetadata = ../../bin/getEnsemblMetadata
checkTransMapSrcDb = ../../bin/checkTransMapSrcDb

diff = diff

test: getGenbankAlignsTests getGenbankMetadataTests getGenbankSeqsTests \
	getEnsemblAlignsTests getEnsemblSeqsTests getEnsemblMetadataTests


###
# getGenbankAligns
###
getGenbankAlignsTests: \
	getGenbankAlignsRnaTest \
	getGenbankAlignsRefSeqTest \
	getGenbankAlignsSplicedEstTest

getGenbankAlignsRnaTest: mkdirs
	${getGenbankAligns} --limit=100 felCat8 mrna output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAligns' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRefs' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

getGenbankAlignsRefSeqTest: mkdirs
	${getGenbankAligns} --limit=100 felCat8 refSeq output/$@.db
	sqlite3 -batch output/$@.db 'select * from srcAligns' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRefs' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

getGenbankAlignsSplicedEstTest: mkdirs
	${getGenbankAligns} --limit=100 felCat8 splicedEst output/$@.db
	sqlite3 -batch output/$@.db 'select * from srcAligns' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRefs' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

###
# getGenbankSeqs
###
getGenbankSeqsTests: \
	getGenbankSeqsRnaTests \
	getGenbankSeqsRefSeqTests

getGenbankSeqsRnaTests: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/rnaAligns.sql
	${getGenbankSeqs} felCat8 mrna output/$@.db
	${checkTransMapSrcDb} felCat8 mrna output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeqs order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids

getGenbankSeqsRefSeqTests: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/rnaAligns.sql
	${getGenbankSeqs} felCat8 refSeq output/$@.db
	${checkTransMapSrcDb} felCat8 refSeq output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq)  as length from srcSeqs order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids

###
# getGenbankMetadata
###
getGenbankMetadataTests: \
	getGenbankMetadataRnaTests \
	getGenbankMetadataRefSeqTests

getGenbankMetadataRnaTests: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/rnaAligns.sql
	${getGenbankMetadata} --testLimit felCat8 mrna output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

getGenbankMetadataRefSeqTests: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/refSeqAligns.sql
	${getGenbankMetadata} --testLimit felCat8 refSeq output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# getEnsemblAligns
##
getEnsemblAlignsTests: \
	getEnsemblAlignsGencodeTest \
	getEnsemblAlignsEnsemblTest

getEnsemblAlignsGencodeTest: mkdirs
	${getEnsemblAligns} --limit=100 hg19 gencode output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAligns' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRefs' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

getEnsemblAlignsEnsemblTest: mkdirs
	${getEnsemblAligns} --limit=100 canFam3 ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header off' 'select * from srcAligns' >output/$@.psl+bin
	sqlite3 -batch output/$@.db 'select * from srcXRefs' >output/$@.xref
	${diff} expected/$@.psl+bin output/$@.psl+bin
	${diff} expected/$@.xref output/$@.xref

##
# getEnsemblSeqs
##
getEnsemblSeqsTests: \
	getEnsemblSeqsGencodeTest \
	getEnsemblSeqsEnsemblTest

getEnsemblSeqsGencodeTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/gencodeAligns.sql
	${getEnsemblSeqs} hg19 gencode output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeqs order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids

getEnsemblSeqsEnsemblTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/gencodeAligns.sql
	${getEnsemblSeqs} hg19 ensembl output/$@.db
	sqlite3 -batch output/$@.db '.header on' 'select name, length(seq) as length from srcSeqs order by name' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids


##
# getEnsemblMetadata
##
getEnsemblMetadataTests: \
	getEnsemblMetadataGencodeTest \
	getEnsemblMetadataEnsemblTest

getEnsemblMetadataGencodeTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/gencodeAligns.sql
	${getEnsemblMetadata} --testLimit hg19 gencode output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta


getEnsemblMetadataEnsemblTest: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/ensemblAligns.sql
	${getEnsemblMetadata} --testLimit canFam3 ensembl output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetadata' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta



##
# Run this to update test data from basic test re
##
update:
	${MAKE} getGenbankAlignsTests diff=true
	sqlite3 -batch output/getGenbankAlignsRefSeqTest.db .dump > input/refSeqAligns.sql
	sqlite3 -batch output/getGenbankAlignsRnaTest.db .dump > input/rnaAligns.sql
	sqlite3 -batch output/getGenbankAlignsSplicedEstTest.db .dump > input/splitEstAligns.sql

	${MAKE} getEnsemblAlignsTests diff=true
	sqlite3 -batch output/getEnsemblAlignsGencodeTest.db .dump > input/gencodeAligns.sql
	sqlite3 -batch output/getEnsemblAlignsEnsemblTest.db .dump > input/ensemblAligns.sql


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
