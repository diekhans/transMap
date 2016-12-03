import os
import sys
from transMap.genomeData import AnnotationType
from pycbio.sys.configInPy import evalConfigFunc

# FIXME: need to track chain type and add to output


class TransMapConf(object):
    """Configuration class for transmap. config file should define
    a function:
       getConfig(<init_keywords_args>)
    The various options are needed only for certain operations and commands; it's an error
    need and not supplied.

    An instance of this class is create by the configuration file, which
    takes set values explicitly and from the command line.
    """
    def __init__(self, configPyFile, dataDir=None, srcHgDb=None, destHgDb=None,
                 annotationType=None, chainType=None, buildTmpDir=None):
        self.configPyFile = configPyFile
        self.dataDir = dataDir
        self.srcHgDb = srcHgDb
        self.destHgDb = destHgDb
        self.annotationType = annotationType
        self.chainType = chainType
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
    def codeRootDir(self):
        return os.path.normpath(os.path.dirname(os.path.dirname(sys.argv[0])))

    @property
    def binDir(self):
        return os.path.join(self.codeRootDir, "bin")

    @property
    def genomeDb(self):
        self.__needOptions("dataDir")
        return os.path.join(self.dataDir, "genome.db")

    @property
    def srcDataDir(self):
        self.__needOptions("dataDir", "srcHgDb")
        return os.path.join(self.dataDir, "src", self.srcHgDb)

    @property
    def srcDb(self):
        self.__needOptions("srcHgDb", "annotationType")
        return os.path.join(self.srcDataDir, "{}.{}.src.db".format(self.srcHgDb, self.annotationType))

    @property
    def mappedDataDir(self):
        self.__needOptions("dataDir", "destHgDb", "srcHgDb")
        return os.path.join(self.dataDir, "mapped", self.destHgDb, self.srcHgDb)

    @property
    def mappedBigPslFile(self):
        self.__needOptions("destHgDb", "srcHgDb", "annotationType")
        return os.path.join(self.mappedDataDir,
                            "{}.{}.{}.mapped.bigPsl".format(self.destHgDb, self.srcHgDb, self.annotationType))

    @property
    def mappingChainsDir(self):
        self.__needOptions("buildTmpDir", "srcHgDb", "destHgDb")
        return os.path.join(self.buildTmpDir, self.destHgDb, self.srcHgDb)

    @property
    def mappingChains(self):
        self.__needOptions("srcHgDb", "destHgDb", "chainType")
        chainFile = "{}.{}.{}.chain".format(self.srcHgDb, self.destHgDb, self.chainType)
        return os.path.join(self.mappingChainsDir, chainFile)

    @property
    def mappingChainsIndexDb(self):
        return self.mappingChains() + ".db"

    @property
    def batchDir(self):
        self.__needOptions("buildTmpDir", "destHgDb", "srcHgDb", "annotationType")
        return os.path.join(self.buildTmpDir, self.destHgDb, self.srcHgDb, str(self.annotationType))

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

    @property
    def jobPreBigPslDir(self):
        return os.path.join(self.batchResutsDir, "parts")

    def getJobPreBigPsl(self, startOid, endOid):
        return os.path.join(self.jobPreBigPslDir, "{}.{}.preBigPsl".format(startOid, endOid))


def transMapConfLoad(configPyFile, dataDir=None, srcHgDb=None, destHgDb=None,
                     annotationType=None, chainType=None, buildTmpDir=None):
    """
    Loads the configuration file which must have define a function getConfig,
    which takes the same signature as TransMapConf.__init__ and returns an instance
    of that class.
    """
    getFuncKwargs = {"dataDir": dataDir,
                     "srcHgDb": srcHgDb,
                     "destHgDb": destHgDb,
                     "annotationType": annotationType,
                     "chainType": chainType,
                     "buildTmpDir": buildTmpDir}
    return evalConfigFunc(configPyFile, getFuncKwargs=getFuncKwargs)
