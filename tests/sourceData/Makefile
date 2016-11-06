
getGenbankAligns = ../../bin/getGenbankAligns
getGenbankSeqs = ../../bin/getGenbankSeqs
getGenbankMetadata = ../../bin/getGenbankMetadata

diff = diff

test: getGenbankAlignsTests getGenbankMetadataTests getGenbankSeqsTests


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
	${diff} expected/$@.psl+bin output/$@.psl+bin

getGenbankAlignsRefSeqTest: mkdirs
	${getGenbankAligns} --limit=100 felCat8 refSeq output/$@.db
	sqlite3 -batch output/$@.db 'select * from srcAligns' >output/$@.psl+bin
	${diff} expected/$@.psl+bin output/$@.psl+bin

getGenbankAlignsSplicedEstTest: mkdirs
	${getGenbankAligns} --limit=100 felCat8 splicedEst output/$@.db
	sqlite3 -batch output/$@.db 'select * from srcAligns' >output/$@.psl+bin
	${diff} expected/$@.psl+bin output/$@.psl+bin

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
	sqlite3 -batch output/$@.db 'select name, length(seq) from srcSeqs' >output/$@.fa.ids
	${diff} expected/$@.fa.ids output/$@.fa.ids

getGenbankSeqsRefSeqTests: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/rnaAligns.sql
	${getGenbankSeqs} felCat8 refSeq output/$@.db
	sqlite3 -batch output/$@.db 'select name, length(seq) from srcSeqs' >output/$@.fa.ids
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
	${getGenbankMetadata} felCat8 mrna output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetaData' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

getGenbankMetadataRefSeqTests: mkdirs
	rm -f output/$@.db
	sqlite3 -batch output/$@.db < input/refSeqAligns.sql
	${getGenbankMetadata} felCat8 refSeq output/$@.db
	sqlite3 -batch output/$@.db  'select * from srcMetaData' >output/$@.meta
	${diff} expected/$@.meta output/$@.meta

##
# Run this to update test data from basic test re
##
update:
	${MAKE} getGenbankAlignsTests diff=true
	sqlite3 -batch output/getGenbankAlignsRefSeqTest.db .dump > input/refSeqAligns.sql
	sqlite3 -batch output/getGenbankAlignsRnaTest.db .dump > input/rnaAligns.sql
	sqlite3 -batch output/getGenbankAlignsSplicedEstTest.db .dump > input/splitEstAligns.sql


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
