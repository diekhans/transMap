import os


class TransMapConf(object):
    """Configuration class for transmap. config file should define
    a function:
       getConfig(<init_keywords_args>)
    The various options are needed only for certain operations and commands; it's an error
    need and not supplied
    """
    def __init__(self, srcDataDir=None, srcHgDb=None, destHgDb=None,
                 annotationType=None, builtTmpDir=None):
        self.srcDataDir = srcDataDir
        self.srcHgDb = srcHgDb
        self.destHgDb = destHgDb
        self.annotationType = annotationType
        self.builtTmpDir = builtTmpDir

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
    def batchDir(self):
        self.__needOptions("buildTmpDir", "srcHgDb", "destHgDb", "annotationType")
        return os.path.join(self.buildTmpDir, self.srcHgDb, self.destHgDb, self.annotationType)

    @property
    def batchParaDir(self):
        return os.path.join(self.batchDir(), "para")

    @property
    def batchParaFile(self):
        return os.path.join(self.batchParaDir(), "transmap.batch")

    @property
    def mappingChains(self):
        chainFile = "{}.{}.{}.chain".format(self.srcHgDb, self.destHgDb, self.annotationType)
        return os.path.join(self.batchDir(), chainFile)

    @property
    def mappingChainsIndexDb(self):
        return self.mappingChains() + ".db"

    @property
    def batchResutsDir(self):
        return os.path.join(self.batchDir(), "results")

    def getJobPreBigBed(self, startOid, endOid):
        return os.path.join(self.batchResutsDir, "{}.{}.preBigBed".format(startOid, endOid))
