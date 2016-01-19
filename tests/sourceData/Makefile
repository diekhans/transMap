
getGenbankAligns = ../../bin/getGenbankAligns
getGenbankSeqs = ../../bin/getGenbankSeqs
getGenbankMetadata = ../../bin/getGenbankMetadata

test: getGenbankAlignsTests getGenbankSeqsTests getGenbankMetadataTests


###
# getGenbankAligns
###
getGenbankAlignsTests: \
	getGenbankAlignsRnaTest \
	getGenbankAlignsRefSeqTest \
	getGenbankAlignsSplicedEstTest

getGenbankAlignsRnaTest: mkdirs
	${getGenbankAligns} felCat8 mrna output/$@.psl output/$@.acc
	diff expected/$@.psl output/$@.psl
	diff expected/$@.acc output/$@.acc

getGenbankAlignsRefSeqTest: mkdirs
	${getGenbankAligns} felCat8 refSeq output/$@.bigPsl output/$@.acc
	bigPslToPsl output/$@.bigPsl stdout | sort >output/$@.psl
	diff expected/$@.psl output/$@.psl
	diff expected/$@.acc output/$@.acc

getGenbankAlignsSplicedEstTest: mkdirs
	${getGenbankAligns} felCat8 splicedEst output/$@.psl.bz2 output/$@.acc.bz2
	bzcat output/$@.psl.bz2 >output/$@.psl
	bzcat output/$@.acc.bz2 >output/$@.acc
	diff expected/$@.psl output/$@.psl
	diff expected/$@.acc output/$@.acc

###
# getGenbankSeqs
###
getGenbankSeqsTests: \
	getGenbankSeqsRnaTests \
	getGenbankSeqsRefSeqTests

getGenbankSeqsRnaTests: mkdirs
	${getGenbankSeqs} felCat8 mrna expected/getGenbankAlignsRnaTest.acc output/$@.fa
	awk '/^>/' output/$@.fa  > output/$@.fa.ids
	diff expected/$@.fa.ids output/$@.fa.ids
	test -e output/$@.fa.fai

getGenbankSeqsRefSeqTests: mkdirs
	${getGenbankSeqs} felCat8 mrna output/getGenbankAlignsRefSeqTest.acc output/$@.fa
	awk '/^>/' output/$@.fa  > output/$@.fa.ids
	diff expected/$@.fa.ids output/$@.fa.ids
	test -e output/$@.fa.fai

###
# getGenbankMetadata
###
getGenbankMetadataTests: \
	getGenbankMetadataRnaTests \
	getGenbankMetadataRefSeqTests

getGenbankMetadataRnaTests: mkdirs
	${getGenbankMetadata} felCat8 mrna expected/getGenbankAlignsRnaTest.acc output/$@.meta
	diff expected/$@.meta output/$@.meta

getGenbankMetadataRefSeqTests: mkdirs
	${getGenbankMetadata} felCat8 mrna expected/getGenbankAlignsRefSeqTest.acc output/$@.meta
	diff expected/$@.meta output/$@.meta

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
