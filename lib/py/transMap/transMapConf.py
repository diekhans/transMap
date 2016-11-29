import os
from transMap.genomeData import AnnotationType

# FIXME: need to track chain type and add to output


class TransMapConf(object):
    """Configuration class for transmap. config file should define
    a function:
       getConfig(<init_keywords_args>)
    The various options are needed only for certain operations and commands; it's an error
    need and not supplied
    """
    def __init__(self, dataDir=None, srcHgDb=None, destHgDb=None,
                 annotationType=None, buildTmpDir=None):
        self.dataDir = dataDir
        self.srcHgDb = srcHgDb
        self.destHgDb = destHgDb
        self.annotationType = annotationType
        self.buildTmpDir = buildTmpDir

        # arguments to pslCDnaFilter
        self.mappedMinQSize = 20
        self.mappedMinCover = 0.20
        self.mappedMaxAligns = 20
        self.mappedGlobalNearBest = 0.01

    def __needOptions(self, *fields):
        for field in fields:
            if getattr(self, field) is None:
                raise Exception("TransMapConf needs {} for this method and it was not specified".format(field))

    @property
    def genomeDataDb(self):
        self.__needOptions("dataDir")
        return os.path.join(self.dataDir, "genome.db")

    @property
    def srcDataDir(self):
        self.__needOptions("dataDir", "srcHgDb")
        return os.path.join(self.dataDir, self.srcHgDb)

    @property
    def srcDb(self):
        self.__needOptions("srcHgDb", "annotationType")
        return os.path.join(self.srcDataDir, "{}.{}.src.db".format(self.srcHgDb, self.annotationType))

    @property
    def mappingChainsDir(self):
        self.__needOptions("buildTmpDir", "srcHgDb", "destHgDb")
        return os.path.join(self.buildTmpDir, self.destHgDb, self.srcHgDb)

    @property
    def mappingChains(self):
        chainFile = "{}.{}.{}.chain".format(self.srcHgDb, self.destHgDb, self.annotationType)
        return os.path.join(self.mappingChainsDir, chainFile)

    @property
    def mappingChainsIndexDb(self):
        return self.mappingChains() + ".db"

    @property
    def batchDir(self):
        self.__needOptions("buildTmpDir", "srcHgDb", "destHgDb", "annotationType")
        return os.path.join(self.buildTmpDir, self.srcHgDb, self.destHgDb, str(self.annotationType))

    @property
    def batchParaDir(self):
        return os.path.join(self.batchDir, "para")

    @property
    def batchParaFile(self):
        return os.path.join(self.batchParaDir, "transmap.batch")

    @property
    def batchResutsDir(self):
        return os.path.join(self.batchDir, "results")

    @property
    def numSeqsPerJob(self):
        "get the number of sequence per batch mapp9ing job for this conf"
        self.__needOptions("srcHgDb", "annotationType")
        if self.annotationType in (AnnotationType.rna, AnnotationType.est):
            return 2500
        else:
            return 500

    def getJobPreBigBed(self, startOid, endOid):
        return os.path.join(self.batchResutsDir, "{}.{}.preBigBed".format(startOid, endOid))
