ROOT = .
.SECONDARY:

# uses in backup
name = transMap-design
docs += transMap-design
#useBibtex = yes
drawings = # output/yyy.pdf

docsPdf = ${docs:%=output/%.pdf}
docsFontchk = ${docs:%=output/%.fontchk}
images = $(wildcard images/*.pdf) $(wildcard images/*.jpg) $(wildcard images/*.png)
depends = ${drawings} ${images}

all: ${docsPdf} ${docsFontchk} ${depends}

ltxopts = -file-line-error-style -output-directory=output
ltxmode =\\nonstopmode\\input

output/%.pdf: %.ltx ${depends}
	@mkdir -p $(dir $@)
	pdflatex ${ltxopts} ${ltxmode} $< </dev/null || (rm -f $@; false)
ifneq (${useBibtex},)
	bibtex output/$*.aux
	pdflatex ${ltxopts} ${ltxmode} $< </dev/null || (rm -f $@; false)
endif
	pdflatex ${ltxopts} ${ltxmode} $< </dev/null || (rm -f $@; false)

output/%.pdf: drawings/%.svg
	@mkdir -p $(dir $@)
	inkscape --export-pdf=$@.tmp $<
	mv -f $@.tmp $@

output/%.fontchk: output/%.pdf
	pdfFontChk $<
	touch $@

output/%.pdf: drawings/%.svg
	@mkdir -p $(dir $@)
	inkscape --export-pdf=$@.tmp $<
	mv -f $@.tmp $@

output/%.pdf: output/%.eps
	epspdf $< $@.tmp.pdf
	mv -f $@.tmp.pdf $@

output/%.eps: drawings/%.gnuplot
	@mkdir -p $(dir $@)
	gnuplot $< >$@.tmp
	mv -f $@.tmp $@

clean:
	rm -rf output _region_.tex

savebak:
	savebak -r hgwdev.cse.ucsc.edu ${name} *.ltx drawings images Makefile

