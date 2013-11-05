ROOT = ../..

build=vertebrate.2012-06-22
paraVer=1
getDataThreads=8

paraTmpDir = para${paraVer}
alignParaDir = ${paraTmpDir}/aln
filterParaDir = ${paraTmpDir}/filt
paraPartDir = ${paraTmpDir}/part

hgLibDir = ${HOME}/compbio/browser/dev/kent/src/hg/lib

host= $(shell hostname)
pid = $(shell echo $$PPID)
mach=$(shell uname -m)
tmpExt = ${host}.${pid}.tmp
tmpExtBz = ${host}.${pid}.tmp.bz2
binDir = ${ROOT}/bin
binDirExec = ${ROOT}/output/bin/${mach}/opt
kentBinDir = ${ROOT}/extern/bin/${mach}
loadedDir = loaded
# this is really under /hive/data/inside/transMap/ aka /cluster/data/transMap/ 
gbdbDir = /gbdb/hgFixed/transMap/${build}

SHELL = bash
export SHELLOPTS=pipefail
ifneq ($(wildcard /data/tmp),)
    export TMPDIR=/data/tmp
else
    export TMPDIR=/scratch/tmp
endif
.SECONDARY:
export PATH:=${binDir}:${binDirExec}:${PATH}

genomeDefs = genomeDefs.pickle
genomeInfo = genomeDefs.info
genomeDists = genomeDefs.dists
mappingsInfo = mappings.info

clades = mammal vertebrate
# Keep public dbs that were transmapped in previous builds
requiredDestDbs = \
    anoCar1 bosTau4 calJac1 canFam2 cavPor3 danRer5 danRer6 dasNov1 echTel1 \
    equCab2 felCat3 fr2 galGal3 gasAcu1 hg18 hg19 loxAfr3 mm9 monDom4 monDom5 \
    ornAna1 oryCun2 oryLat2 otoGar1 panTro2 petMar1 ponAbe2 rheMac2 rn4 taeGut1 \
    tetNig1 tetNig2 xenTro2

synOpts=--synDist=0.5 --useAvailableRBest
cladeOpts = $(foreach c,${clades},--includeClade=${c})
destDbOpts =   $(foreach d,${requiredDestDbs},--otherDestDb=${d})
dbOpts= ${cladeOpts} ${destDbOpts}

# rn5 is deadend
# dasNov2 is stalled
excludeDbs = venter1 canFamPoodle1 rn5.obsolete dasNov2
excludeDbOpts = $(foreach d,${excludeDbs},--excludeDb=${d})
trees =	--tree=/hive/data/mm9/bed/multiz30way/mm9.30way.mod \
	--tree=/hive/data/mm10/bed/multiz60way/mm10.60way.nh \
	--tree=/hive/data/hg19/bed/multiz46way/46way.nh  \
	--tree=/hive/data/hg18/bed/multiz44way/44way.nh \
	--tree=/hive/data/petMar1/bed/multiz6way/downloads/phastCons6way/6way.mod
defsOpts= ${cladeOpts} ${excludeDbOpts} ${trees}


remConf=${HOME}/.hg.rem.conf
localConf=${HOME}/.hg.conf

ifeq (${HGDB_CONF},)
    export HGDB_CONF=${remConf}
endif

all:

getData:  ${genomeDefs} ${mappingsInfo}

defs: ${genomeDefs}


####
${genomeDefs}:
	${binDir}/genomeDefsMk ${defsOpts} --info=${genomeInfo} --chainDists=${genomeDists} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${mappingsInfo}: ${genomeDefs}
	${binDir}/getData ${dbOpts} --threads=${getDataThreads} ${genomeDefs} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

batchSetup:
	@mkdir -p  ${alignParaDir} ${filterParaDir}  ${paraPartDir}
	${binDir}/mkTransMapJobs ${dbOpts} ${synOpts} ${genomeDefs} ${paraDir} ${alignParaDir} ${filterParaDir} ${paraPartDir} aligns

cleanDefs:
	rm -f ${genomeDefs} ${genomeInfo} ${genomeDists}

batchCleanup:
	rm -rf ${alignParaDir} ${filterParaDir}  ${paraPartDir}

destDbs = $(notdir $(wildcard aligns/*))
srcTypes = refSeq mrna splicedEst ucscGenes 
srcTypeTblPartMap_refSeq = RefSeq
srcTypeTblPartMap_mrna = MRna
srcTypeTblPartMap_splicedEst = SplicedEst
srcTypeTblPartMap_ucscGenes = UcscGenes


###
# database loading
###
loadDb: loadTracks loadMeta
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
	bzcat $(wildcard aligns/${destDb}/*/$*.psl.bz2) | sort -k 14,14 -k 16,16n  | hgLoadPsl ${destDb} -table=transMapAln${srcTypeTblPartMap_$*} stdin
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
loadSrcMeta: ${loadedDir}/load.meta.src.${srcType}.done 

${loadedDir}/load.meta.src.${srcType}.done:
	@mkdir -p $(dir $@)
	mkTransMapSrc $(wildcard data/*/${srcType}.psl.bz2) | hgLoadSqlTab hgFixed transMapSrc${srcTypeTbl} ${hgLibDir}/transMapSrc.sql stdin
	touch $@

# srcType= srcTypeTbl=
loadGeneMeta: ${loadedDir}/load.meta.gene.${srcType}.done

geneAlignedIdTmp = ${srcType}.alignedId.${tmpExt}

${loadedDir}/load.meta.gene.${srcType}.done:
	@mkdir -p $(dir $@)
ifneq (${srcType},splicedEst)
	bzcat $(wildcard data/*/${srcType}.psl.bz2) | cut -f 10 >${geneAlignedIdTmp}
	mkTransMapGene ${geneAlignedIdTmp} $(wildcard data/*/${srcType}.meta.bz2) | hgLoadSqlTab hgFixed transMapGene${srcTypeTbl} ${hgLibDir}/transMapGene.sql stdin
	rm -f ${geneAlignedIdTmp}
endif
	touch $@

# srcType= srcTypeTbl=
loadSeq: ${loadedDir}/load.meta.seq.${srcType}.done

${loadedDir}/load.meta.seq.${srcType}.done: ${gbdbDir}/${srcType}.fa
	@mkdir -p $(dir $@)
	hgsql -e 'drop table if exists transMapSeq${srcTypeTbl}' hgFixed
	hgsql -e 'drop table if exists transMapExtFile${srcTypeTbl}' hgFixed
	hgLoadSeq -seqTbl=transMapSeq${srcTypeTbl} -extFileTbl=transMapExtFile${srcTypeTbl} hgFixed $<
	rm -f transMapSeq${srcTypeTbl}.tab
	touch $@

${gbdbDir}/${srcType}.fa: $(wildcard data/*/${srcType}.fa.bz2)
	@mkdir -p $(dir $@)
	bzcat $^ | faFilter -uniq stdin $@.tmp
	mv -f $@.tmp $@

####
allJoiner=${HOME}/compbio/browser/dev/kent/src/hg/makeDb/schema/all.joiner
genomeJoinerCheckIds = transMapAlnInfoRefSeqId transMapAlnInfoUcscGenesId transMapAlnInfoMRnaId transMapAlnInfoSplicedEstId
hgFixedJoinerCheckIds = transMapSrcRefSeqId transMapSrcUcscGenesId transMapSrcMRnaId transMapSrcSplicedEstId transMapExtFileRefSeqId transMapExtFileUcscGenesId transMapExtFileMRnaId transMapExtFileSplicedEstId  transMapSrcInfoRefSeqId transMapSrcInfoUcscGenesId transMapSrcInfoMRnaId transMapSrcInfoSplicedEstId transMapSrcAlnUcscGenesId transMapSrcAlnRefSeqId transMapSrcAlnSplicedEstId
joinerCheckIds = ${genomeJoinerCheckIds} ${hgFixedJoinerCheckIds}
joinerCheckOutDir = jchecks

# must run for all databases, or cross-database links are not checked
runJoinerCheck: $(joinerCheckIds:%=${joinerCheckOutDir}/%.jcheck) \

# -keys must be specified with -identifier= 
${joinerCheckOutDir}/%.jcheck:
	@mkdir -p $(dir $@)
	joinerCheck -keys -identifier=$* ${allJoiner} >&$@.tmp
	mv -f $@.tmp $@

####
# generate script to pack tables
mkPackScript:
	cp -f misc/transMapPack.head misc/transMapPack
	mkPackScript ${destDbs} hgFixed >>misc/transMapPack
	chmod +x misc/transMapPack

