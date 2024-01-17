
pyprogs = $(shell file -F $$'\t' bin/* tests/*/bin/* | awk '/Python script/{print $$1}')
export PYTHONWARNINGS=always

all:

test:
	(cd tests && ${MAKE} test)

lint:
	python3 -m flake8 --color=never lib/py/transMap ${pyprogs}

clean:
	(cd tests && ${MAKE} clean)
	rm -f lib/py/transMap/*.pyc
