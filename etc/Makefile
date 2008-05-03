ROOT = ../..
pid = $(shell echo $$PPID)
tmpExt = ${HOST}.${pid}.tmp
tmpExtBz = ${HOST}.${pid}.tmp.bz2
binDir = ${ROOT}/bin

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

paraDir = para3
scratchOutDir = /san/sanvol1/scratch/markd/transMap/${para1}
jobSpec = ${paraDir}/jobs.para


export HGDB_CONF=${HOME}/.hg.rem.conf

all: ${genomeDefs} ${mappingsInfo}

defs: ${genomeDefs}

batchSetup: ${jobSpec}

${genomeDefs}:
	${binDir}/genomeDefsMk --info=${genomeInfo} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${mappingsInfo}: ${genomeDefs}
	${binDir}/getData ${clades} ${onlyDestDbs} ${otherDestDbs} ${skipDbs}  ${genomeDefs} $@.${tmpExt}
	mv -f $@.${tmpExt} $@


${jobSpec}:
	@mkdir -p ${paraDir} ${scratchOutDir}
	${binDir}/mkTransMapJobs ${genomeDefs} ${paraDir} ${scratchOutDir}
