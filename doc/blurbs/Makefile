
images = transMap-xspecies transMap-transcriptome

all: ${images:%=%.pdf}

%.pdf: %.svg
	svg2pdf $< $@


clean:
	rm -f ${images:%=%.pdf}
