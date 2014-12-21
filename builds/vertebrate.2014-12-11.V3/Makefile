ROOT = ../..

version = V3
paraVer = 1
getDataThreads = 16

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
gbdbDir = /gbdb/hgFixed/transMap/${version}

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
genomeInfoTsv = genomeDefs.tsv
genomeChainTsv = genomeDefs.chain.tsv
mappingsInfo = mappings.info

# this comes from hgcentraltest, which is different than public
clades = mammal primate vertebrate rodent
# Assemblies to include even if they are not the latest.
requiredDestDbs = mm9 hg19 rn5 felCat5

synOpts=--synDist=1.0 --useAvailableRBest
cladeOpts = $(foreach c,${clades},--includeClade=${c})
destDbOpts =   $(foreach d,${requiredDestDbs},--otherDestDb=${d})
dbOpts= ${cladeOpts} ${destDbOpts}

# rn5.obsolete is dead-end
# dasNov2 is stalled
excludeDbs = venter1 canFamPoodle1 rn5.obsolete dasNov2
excludeDbOpts = $(foreach d,${excludeDbs},--excludeDb=${d})
trees =	--tree=/hive/data/genomes/mm9/bed/multiz30way/mm9.30way.mod \
	--tree=/hive/data/genomes/mm10/bed/multiz60way/mm10.60way.nh \
	--tree=/hive/data/genomes/hg19/bed/multiz46way/46way.nh  \
	--tree=/hive/data/genomes/hg18/bed/multiz44way/44way.nh \
	--tree=/hive/data/genomes/hg19/bed/multiz100way/hg19.100way.nh \
	--tree=/hive/data/genomes/petMar1/bed/multiz6way/downloads/phastCons6way/6way.mod
defsOpts= ${cladeOpts} ${excludeDbOpts} ${trees}


remConf=${HOME}/.hg.rem.conf
localConf=${HOME}/.hg.conf

ifeq (${HGDB_CONF},)
    export HGDB_CONF=${remConf}
endif

all: getDefs getData

getDefs: ${genomeDefs}

getData:  ${genomeDefs} ${mappingsInfo}

####
# geneome definitions
${genomeDefs}:
	${binDir}/genomeDefsMk ${defsOpts} --genomeInfo=${genomeInfoTsv} --chainInfo=${genomeChainTsv} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####
# obtain data
${mappingsInfo}: ${genomeDefs}
	${binDir}/getData ${dbOpts} --threads=${getDataThreads} ${genomeDefs} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####
# mapping

batchSetup:
	@mkdir -p  ${alignParaDir} ${filterParaDir}  ${paraPartDir}
	${binDir}/mkTransMapJobs ${dbOpts} ${synOpts} ${genomeDefs} ${paraDir} ${alignParaDir} ${filterParaDir} ${paraPartDir} aligns


defsClean:
	rm -f ${genomeDefs} ${genomeInfoTsv} ${genomeChainTsv}

batchClean:
	rm -rf ${alignParaDir} ${filterParaDir}  ${paraPartDir}

destDbs = $(notdir $(wildcard aligns/*))
srcTypes = refSeq mrna splicedEst ucscGenes ensembl
srcTypeTblPartMap_refSeq = RefSeq
srcTypeTblPartMap_mrna = MRna
srcTypeTblPartMap_splicedEst = SplicedEst
srcTypeTblPartMap_ucscGenes = UcscGenes
srcTypeTblPartMap_ensembl = Ensembl


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
	bzcat $(wildcard aligns/${destDb}/*/$*.psl.bz2) | sort -k 14,14 -k 16,16n  | hgLoadPsl ${destDb} -table=transMapAln${srcTypeTblPartMap_$*}${version} stdin
	touch $@

${loadedDir}/load.${destDb}.info.%.done:
	@mkdir -p $(dir $@)
	mkTransMapInfo $(wildcard aligns/${destDb}/*/$*.mapinfo.bz2) | hgLoadSqlTab ${destDb} transMapInfo${srcTypeTblPartMap_$*}${version} ${hgLibDir}/transMapInfo.sql stdin
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
	mkTransMapSrc $(wildcard data/*/${srcType}.psl.bz2) | hgLoadSqlTab hgFixed transMapSrc${srcTypeTbl}${version} ${hgLibDir}/transMapSrc.sql stdin
	touch $@

# srcType= srcTypeTbl=
loadGeneMeta: ${loadedDir}/load.meta.gene.${srcType}.done

geneAlignedIdTmp = ${srcType}.alignedId.${tmpExt}

${loadedDir}/load.meta.gene.${srcType}.done:
	@mkdir -p $(dir $@)
ifneq (${srcType},splicedEst)
	bzcat $(wildcard data/*/${srcType}.psl.bz2) | cut -f 10 >${geneAlignedIdTmp}
	mkTransMapGene ${geneAlignedIdTmp} $(wildcard data/*/${srcType}.meta.bz2) | hgLoadSqlTab hgFixed transMapGene${srcTypeTbl}${version} ${hgLibDir}/transMapGene.sql stdin
	rm -f ${geneAlignedIdTmp}
endif
	touch $@

# srcType= srcTypeTbl=
loadSeq: ${loadedDir}/load.meta.seq.${srcType}.done

${loadedDir}/load.meta.seq.${srcType}.done: ${gbdbDir}/${srcType}.fa
	@mkdir -p $(dir $@)
	hgsql -e 'drop table if exists transMapSeq${srcTypeTbl}${version}' hgFixed
	hgsql -e 'drop table if exists transMapExtFile${srcTypeTbl}${version}' hgFixed
	hgLoadSeq -seqTbl=transMapSeq${srcTypeTbl}${version} -extFileTbl=transMapExtFile${srcTypeTbl}${version} hgFixed $<
	rm -f transMapSeq${srcTypeTbl}${version}.tab
	touch $@

${gbdbDir}/${srcType}.fa: $(wildcard data/*/${srcType}.fa.bz2)
	@mkdir -p $(dir $@)
	bzcat $^ | faFilter -uniq stdin $@.tmp
	mv -f $@.tmp $@

####
allJoiner=${HOME}/compbio/browser/dev/kent/src/hg/makeDb/schema/all.joiner
genomeJoinerCheckIds = transMapAlnInfoRefSeq${version}Id transMapAlnInfoUcscGenes${version}Id transMapAlnInfoMRna${version}Id transMapAlnInfoSplicedEst${version}Id
hgFixedJoinerCheckIds = transMapSrcRefSeq${version}Id transMapSrcUcscGenes${version}Id transMapSrcMRna${version}Id transMapSrcSplicedEst${version}Id transMapExtFileRefSeq${version}Id transMapExtFileUcscGenes${version}Id transMapExtFileMRna${version}Id transMapExtFileSplicedEst${version}Id  transMapSrcInfoRefSeq${version}Id transMapSrcInfoUcscGenes${version}Id transMapSrcInfoMRna${version}Id transMapSrcInfoSplicedEst${version}Id transMapSrcAlnUcscGenes${version}Id transMapSrcAlnRefSeq${version}Id transMapSrcAlnSplicedEst${version}Id
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
# generate list of databases to push
mkPushList:


####
# generate script to pack tables
mkPackScript:
	cp -f misc/transMapPack.head misc/transMapPack
	mkPackScript ${version} ${destDbs} hgFixed >>misc/transMapPack
	chmod +x misc/transMapPack

