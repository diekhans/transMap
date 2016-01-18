
getGenbankAligns = ../../bin/getGenbankAligns


test: getGenbankAlignsTests


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
	${getGenbankAligns} felCat8 refSeq output/$@.psl output/$@.acc
	diff expected/$@.psl output/$@.psl
	diff expected/$@.acc output/$@.acc

getGenbankAlignsSplicedEstTest: mkdirs
	${getGenbankAligns} felCat8 splicedEst output/$@.psl output/$@.acc
	diff expected/$@.psl output/$@.psl
	diff expected/$@.acc output/$@.acc

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
