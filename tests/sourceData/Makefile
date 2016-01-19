
getGenbankAligns = ../../bin/getGenbankAligns
getGenbankSeqs = ../../bin/getGenbankSeqs


test: getGenbankAlignsTests getGenbankSeqsTests


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

getGenbankSeqsRnaTests: mkdirs getGenbankAlignsRnaTest
	${getGenbankSeqs} felCat8 mrna output/getGenbankAlignsRnaTest.acc output/$@.fa
	awk '/^>/' output/$@.fa  > output/$@.fa.ids
	diff expected/$@.fa.ids output/$@.fa.ids
	test -e output/$@.fa.fai

getGenbankSeqsRefSeqTests: mkdirs getGenbankAlignsRefSeqTest
	${getGenbankSeqs} felCat8 mrna output/getGenbankAlignsRefSeqTest.acc output/$@.fa
	awk '/^>/' output/$@.fa  > output/$@.fa.ids
	diff expected/$@.fa.ids output/$@.fa.ids
	test -e output/$@.fa.fai

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
