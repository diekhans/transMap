ROOT = ../..

build =vertebrate.2008-05-20
paraVer=2

alignParaDir = paraAln${paraVer}
filterParaDir = paraFilt${paraVer}

hgLibDir = ${HOME}/compbio/genbank/kent/src/hg/lib

pid = $(shell echo $$PPID)
tmpExt = ${HOST}.${pid}.tmp
tmpExtBz = ${HOST}.${pid}.tmp.bz2
binDir = ${ROOT}/bin
loadedDir = loaded
# this is really under /cluster/data/transMap/
gbdbDir = /gbdb/hgFixed/transMap/${build}

SHELL = bash3.2
export SHELLOPTS=pipefail
export TMPDIR=/scratch/tmp
.SECONDARY:
export PATH:=${binDir}:${PATH}

genomeDefs = genomeDefs.pickle
genomeInfo = genomeDefs.info
mappingsInfo = mappings.info

clades = --clade=mammal --clade=vertebrate
#onlyDestDbs = --onlyDestDb=XXX 
otherDestDbs = --otherDestDb=hg17 --otherDestDb=mm8
skipDbs = --skipDb=venter1
dbOpts= ${clades} ${onlyDestDbs} ${otherDestDbs} ${skipDbs}

scartchRootDir = /san/sanvol1/scratch/markd/transMap
scratchOutDir = ${scartchRootDir}/${build}

remConf=${HOME}/.hg.rem.conf
localConf=${HOME}/.hg.conf

ifeq (${HGDB_CONF},)
    export HGDB_CONF=${remConf}
endif

all:

getData:  ${genomeDefs} ${mappingsInfo}

defs: ${genomeDefs}

batchSetup:
	@mkdir -p  ${alignParaDir} ${filterParaDir}  ${scratchOutDir}
	${binDir}/mkTransMapJobs ${dbOpts} ${genomeDefs} ${paraDir} ${alignParaDir} ${filterParaDir} ${scratchOutDir} aligns

${genomeDefs}:
	${binDir}/genomeDefsMk --info=${genomeInfo} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${mappingsInfo}: ${genomeDefs}
	${binDir}/getData ${dbOpts} ${genomeDefs} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

copyAlns:
	@mkdir -p aligns
	rsync -a --whole-file ${scratchOutDir}/ aligns

destDbs = $(notdir $(wildcard aligns/*))
srcTypes = refSeq mrna splicedEst ucscGenes 
srcTypeTblPartMap_refSeq = RefSeq
srcTypeTblPartMap_mrna = MRna
srcTypeTblPartMap_splicedEst = SplicedEst
srcTypeTblPartMap_ucscGenes = UcscGenes


###
# database loading
###
loadTracks: $(destDbs:%=%.trackLoad)
loadMeta: $(srcTypes:%=%.metaLoad)

#####
# per-genome load  %==destDb
###
%.trackLoad:
	@HGDB_CONF=${localConf} ${MAKE} loadTrackTypes destDb=$*

# params: destDb=
destDbTypes = $(filter ${srcTypes},$(sort $(basename $(basename $(notdir $(wildcard aligns/${destDb}/*/*.psl.bz2))))))

loadTrackTypes: \
	$(destDbTypes:%=${loadedDir}/load.${destDb}.aln.%.done) \
	$(destDbTypes:%=${loadedDir}/load.${destDb}.info.%.done)

${loadedDir}/load.${destDb}.aln.%.done:
	@mkdir -p $(dir $@)
	bzcat $(wildcard aligns/${destDb}/*/$*.psl.bz2) | hgLoadPsl ${destDb} -table=transMapAln${srcTypeTblPartMap_$*} stdin
	touch $@

${loadedDir}/load.${destDb}.info.%.done:
	@mkdir -p $(dir $@)
	mkTransMapInfo $(wildcard aligns/${destDb}/*/$*.mapinfo.bz2) | hgLoadSqlTab ${destDb} transMapInfo${srcTypeTblPartMap_$*} ${hgLibDir}/transMapInfo.sql stdin
	touch $@

#####
# meta-data load  %=srcType
###
%.metaLoad:
	@HGDB_CONF=${localConf} ${MAKE} loadSrcMeta srcType=$* srcTypeTbl=${srcTypeTblPartMap_$*}
	@HGDB_CONF=${localConf} ${MAKE} loadGeneMeta srcType=$* srcTypeTbl=${srcTypeTblPartMap_$*}
	@HGDB_CONF=${localConf} ${MAKE} loadSeq srcType=$* srcTypeTbl=${srcTypeTblPartMap_$*}

# srcType= srcTypeTbl=
loadSrcMeta: ${loadedDir}/load.src.${srcType}.done 

${loadedDir}/load.src.${srcType}.done:
	@mkdir -p $(dir $@)
	mkTransMapSrc $(wildcard data/*/${srcType}.psl.bz2) | hgLoadSqlTab hgFixed transMapSrc${srcTypeTbl} ${hgLibDir}/transMapSrc.sql stdin
	touch $@

# srcType= srcTypeTbl=
loadGeneMeta: ${loadedDir}/load.gene.${srcType}.done

geneAlignedIdTmp = ${srcType}.alignedId.${tmpExt}

${loadedDir}/load.gene.${srcType}.done:
	@mkdir -p $(dir $@)
ifneq (${srcType},splicedEst)
	bzcat $(wildcard data/*/${srcType}.psl.bz2) | cut -f 10 >${geneAlignedIdTmp}
	mkTransMapGene ${geneAlignedIdTmp} $(wildcard data/*/${srcType}.meta.bz2) | hgLoadSqlTab hgFixed transMapGene${srcTypeTbl} ${hgLibDir}/transMapGene.sql stdin
	rm -f ${geneAlignedIdTmp}
endif
	touch $@

# srcType= srcTypeTbl=
loadSeq: ${loadedDir}/load.seq.${srcType}.done

${loadedDir}/load.seq.${srcType}.done: ${gbdbDir}/${srcType}.fa
	@mkdir -p $(dir $@)
	hgsql -e 'drop table if exists transMapSeq${srcTypeTbl}' hgFixed
	hgsql -e 'drop table if exists transMapExtFile${srcTypeTbl}' hgFixed
	hgLoadSeq -seqTbl=transMapSeq${srcTypeTbl} -extFileTbl=transMapExtFile${srcTypeTbl} hgFixed $<
	touch $@

${gbdbDir}/${srcType}.fa: $(wildcard data/*/${srcType}.fa.bz2)
	@mkdir -p $(dir $@)
	bzcat $^ | faFilter -uniq stdin $@.tmp
	mv -f $@.tmp $@
