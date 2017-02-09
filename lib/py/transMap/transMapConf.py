import os
import sys
from transMap.genomeData import AnnotationType
from pycbio.sys.configInPy import evalConfigFunc
from transMap.genbankConf import GenbankConf
from transMap.genomeData import ChainType

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
                 annotationType=None, chainType=None, buildTmpDir=None,
                 version="", batchGen=1):
        self.configPyFile = configPyFile
        self.dataDir = dataDir
        self.srcHgDb = srcHgDb
        self.destHgDb = destHgDb
        self.annotationType = annotationType
        self.chainType = chainType
        self.buildTmpDir = buildTmpDir
        self.version = version
        self.batchGen = batchGen

        # ucsc browser data
        self.hgCentralDb = "hgcentraltest"

        # genbank
        self.genbankConfRa = GenbankConf.stdConfRaFile

        # install
        self.gbdbDir = "/gbdb"
        
        # selecting genomes
        self.clades = frozenset(["vertebrate", "mammal"])

        # required older versions of some genomes
        self.requiredPreviousDestHgDbs = frozenset([
            "hg19"])

        # types of chains for mappings between clades each entry is a set of
        # clades and an ordered tuple of preferred chain types
        self.cladePreferedChains = {
            frozenset(["mammal", "mammal"]): (ChainType.syn, ChainType.rbest, ChainType.all),
            frozenset(["vertebrate", "mammal"]): (ChainType.rbest, ChainType.all),
            frozenset(["vertebrate", "vertebrate"]): (ChainType.rbest, ChainType.all),
        }

        # arguments to pslCDnaFilter
        self.mappedMinQSize = 20
        self.mappedMinCover = 0.20
        self.mappedMaxAligns = 20
        self.mappedGlobalNearBest = 0.01

        self.paraHost = "ku"

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

    def getMappedDataDirForDestHgDb(self, destHgDb):
        self.__needOptions("dataDir")
        return os.path.join(self.dataDir, "mapped", destHgDb)

    @property
    def mappedDataDir(self):
        self.__needOptions("dataDir", "destHgDb")
        return self.getMappedDataDirForDestHgDb(self.destHgDb)

    def __bigPslFileBasename(self, destHgDb, annotationType):
        return "{}.{}.transMap{}.bigPsl".format(destHgDb, annotationType, self.version)
            
    def getMappedBigPslFileForDestHgDb(self, destHgDb, annotationType):
        return os.path.join(self.getMappedDataDirForDestHgDb(destHgDb),
                            self.__bigPslFileBasename(destHgDb, annotationType))

    def getMappedGbdbDir(self, destHgDb):
        return os.path.join(self.gbdbDir, destHgDb, "transMap", self.version)
    
    def getMappedGbdbBigPslFile(self, destHgDb, annotationType):
        return os.path.join(self.getMappedGbdbDir(destHgDb),
                            self.__bigPslFileBasename(destHgDb, annotationType))

    @property
    def mappedBigPslFile(self):
        self.__needOptions("destHgDb", "annotationType")
        return self.getMappedBigPslFileForDestHgDb(self.destHgDb, self.annotationType)

    @property
    def mappingChainsDir(self):
        self.__needOptions("buildTmpDir", "srcHgDb", "destHgDb")
        return os.path.join(self.buildTmpDir, "chains", self.destHgDb, self.srcHgDb)

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
        return os.path.join(self.buildTmpDir, "batches", self.destHgDb, self.srcHgDb, str(self.annotationType))

    @property
    def batchParaDir(self):
        return os.path.join(self.batchDir, "para{}".format(self.batchGen))

    @property
    def batchParaFile(self):
        return os.path.join(self.batchParaDir, "transmap.batch")

    @property
    def batchResutsDir(self):
        return os.path.join(self.batchDir, "results")

    @property
    def numSeqsPerJob(self):
        "get the number of sequence per batch mapping job for this conf"
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

    def getBatchSrcHgDbPreBigPsl(self, srcHgDb):
        self.__needOptions("buildTmpDir", "destHgDb", "annotationType")
        return os.path.join(self.buildTmpDir, "results", self.destHgDb,
                            "{}.{}.{}.preBigPsl".format(self.destHgDb, srcHgDb, self.annotationType))

    @property
    def batchPreBigPsl(self):
        self.__needOptions("buildTmpDir", "srcHgDb", "destHgDb", "annotationType")
        return self.getBatchSrcHgDbPreBigPsl(self.srcHgDb)

    @property
    def bigTransMapAsPath(self):
        return os.path.expanduser("~/kent/src/hg/lib/bigTransMap.as")

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
