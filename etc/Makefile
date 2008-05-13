ROOT = ../..

paraDir = para9

hgLibDir = ${HOME}/compbio/genbank/kent/src/hg/lib

pid = $(shell echo $$PPID)
tmpExt = ${HOST}.${pid}.tmp
tmpExtBz = ${HOST}.${pid}.tmp.bz2
binDir = ${ROOT}/bin
loadedDir = loaded

SHELL = bash3.2
export SHELLOPTS=pipefail
export TMPDIR=/scratch/tmp
.SECONDARY:
export PATH:=${binDir}:${PATH}

genomeDefs = genomeDefs.pickle
genomeInfo = genomeDefs.info
mappingsInfo = mappings.info

clades = --clade=mammal --clade=vertebrate
onlyDestDbs = 
otherDestDbs = --otherDestDb=hg17 --otherDestDb=mm8
skipDbs = --skipDb=venter1
dbOpts= ${clades} ${onlyDestDbs} ${otherDestDbs} ${skipDbs}

scratchOutDir = /san/sanvol1/scratch/markd/transMap
jobSpec = ${paraDir}/jobs.para

remConf=${HOME}/.hg.rem.conf
localConf=${HOME}/.hg.conf

ifeq (${HGDB_CONF},)
    export HGDB_CONF=${remConf}
endif

all: ${genomeDefs} ${mappingsInfo}

defs: ${genomeDefs}

batchSetup: ${jobSpec}

${genomeDefs}:
	${binDir}/genomeDefsMk --info=${genomeInfo} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${mappingsInfo}: ${genomeDefs}
	${binDir}/getData ${dbOpts} ${genomeDefs} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${jobSpec}:
	@mkdir -p ${paraDir} ${scratchOutDir}
	${binDir}/mkTransMapJobs ${dbOpts} ${genomeDefs} ${paraDir} ${scratchOutDir}

copyAlns:
	@mkdir -p aligns
	rsync -a --whole-file ${scratchOutDir}/ aligns

destDbs = $(notdir $(wildcard aligns/*))
srcTypes = refSeq
# FIXME: mrna splicedEst
# FIXME ucscGenes
srcTypeTblPartMap_refSeq = RefSeq
srcTypeTblPartMap_mrna = MRna
srcTypeTblPartMap_splicedEst = SplicedEst
srcTypeTblPartMap_ucscGenes = UcscGenes


loadTracks: $(destDbs:%=%.trackLoad) $(srcTypes:%=%.metaLoad)

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
ifneq (${srcType}, splicedEst)
	@HGDB_CONF=${localConf} ${MAKE} loadGeneMeta srcType=$* srcTypeTbl=${srcTypeTblPartMap_$*}
endif

# srcType= srcTypeTbl=
loadSrcMeta: ${loadedDir}/load.src.${srcType}.done

${loadedDir}/load.src.${srcType}.done:
	@mkdir -p $(dir $@)
	mkTransMapSrc $(wildcard data/*/${srcType}.psl.bz2) | hgLoadSqlTab hgFixed transMapSrc${srcTypeTbl} ${hgLibDir}/transMapSrc.sql stdin
	touch $@

# srcType= srcTypeTbl=
loadGeneMeta: ${loadedDir}/load.gene.${srcType}.done

${loadedDir}/load.gene.${srcType}.done:
	@mkdir -p $(dir $@)
	mkTransMapGene $(wildcard data/*/${srcType}.meta.bz2) | hgLoadSqlTab hgFixed transMapGene{srcTypeTbl} ${hgLibDir}/transMapGene.sql stdin
	touch $@
