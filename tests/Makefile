all:

test: sourceData.test

%.test:
	(cd $* && ${MAKE} test)

clean: sourceData.clean

%.clean:
	(cd $* && ${MAKE} clean)
