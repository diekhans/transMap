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
    - srcTwoBitPath - a list of path names to search to find srcDb twobit, with the
      database formatted as {db}. A path allows for checking on cluster local
      storage.  A single string may also be supplied.
    - srcChromSizes - location of srcDb chrom sizes file.
    - destTwoBitPathPat - a list of path names to search to find destDb twobit, with the
      database formatted as {db}. A path allows for checking on cluster local
      storage.  A single string may also be supplied.
    - destChromSizesPat - location of destDb chrom sizes file.

    """
    def __init__(self, configPyFile, *, dataRootDir=None, srcHgDb=None, destHgDb=None,
                 annotationType=None, chainType=None, version=None, batchGen=None,
                 prevVersion=None, prevDataRootDir=None,
                 srcTwoBitPathPat=None, srcChromSizesPat=None,
                 destTwoBitPathPat=None, destChromSizesPat=None, numSeqsPerJob=None,
                 paraHost="ku"):
        self.configPyFile = configPyFile
        self.dataRootDir = dataRootDir
        self.srcHgDb = srcHgDb
        self.destHgDb = destHgDb
        self.annotationType = annotationType
        self.chainType = chainType
        self.version = version
        self.batchGen = batchGen
        self.prevVersion = prevVersion
        self.prevDataRootDir = prevDataRootDir
        self.srcTwoBitPathPat = srcTwoBitPathPat
        self.srcChromSizesPat = srcChromSizesPat
        self.destTwoBitPathPat = destTwoBitPathPat
        self.destChromSizesPat = destChromSizesPat
        self.fixedNumSeqsPerJob = numSeqsPerJob

        # ucsc browser data
        self.hgCentralDb = "hgcentraltest"

        # genome data defaults
        _defaultTwoBitPathPat = ("/scratch/data/{db}/{db}.2bit",
                                 "/hive/data/genomes/{db}/{db}.2bit")
        _defaultChromSizesPat = "/hive/data/genomes/{db}/chrom.sizes"
        if self.srcTwoBitPathPat is None:
            self.srcTwoBitPathPat = _defaultTwoBitPathPat
        elif isinstance(self.srcTwoBitPathPat, str):
            self.srcTwoBitPathPat = [self.srcTwoBitPathPat]
        if self.srcChromSizesPat is None:
            self.srcChromSizesPat = _defaultChromSizesPat
        if self.destTwoBitPathPat is None:
            self.destTwoBitPathPat = _defaultTwoBitPathPat
        elif isinstance(self.destTwoBitPathPat, str):
            self.destTwoBitPathPat = [self.destTwoBitPathPat]
        if self.destChromSizesPat is None:
            self.destChromSizesPat = _defaultChromSizesPat

        # genbank
        self.genbankConfRa = GenbankConf.stdConfRaFile

        # install
        self.gbdbDir = "/gbdb"

        # selecting genomes
        self.clades = frozenset(["vertebrate", "mammal"])

        # required older versions of some genomes.  This is load in a lacy manned and obtained
        # from prevDataRootDir
        self.requiredPreviousDestHgDbsCache = None

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

        self.paraHost = paraHost

    def _needOptions(self, *fields):
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
    def dataDir(self):
        self._needOptions("dataRootDir")
        return os.path.normpath(os.path.join(self.dataRootDir, "data"))

    @property
    def buildDir(self):
        self._needOptions("dataRootDir")
        return os.path.join(self.dataRootDir, "build")

    @property
    def genomeDb(self):
        return os.path.join(self.dataDir, "genome.db")

    def getSrcDataDir(self, srcHgDb):
        return os.path.join(self.dataDir, "src", srcHgDb)

    @property
    def srcDataDir(self):
        self._needOptions("srcHgDb")
        return self.getSrcDataDir(self.srcHgDb)

    def getSrcDb(self, srcHgDb, annotationType):
        return os.path.join(self.getSrcDataDir(srcHgDb), "{}.{}.src.db".format(srcHgDb, annotationType))

    @property
    def srcDb(self):
        self._needOptions("srcHgDb", "annotationType")
        return self.getSrcDb(self.srcHgDb, self.annotationType)

    @property
    def resultsDir(self):
        self._needOptions("dataRootDir")
        return os.path.normpath(os.path.join(self.dataRootDir, "results"))

    def getMappedResultsForDestHgDb(self, destHgDb):
        return os.path.join(self.resultsDir, "mapped", destHgDb)

    @property
    def mappedResultsDir(self):
        self._needOptions("destHgDb")
        return self.getMappedResultsForDestHgDb(self.destHgDb)

    def _bigPslFileBasename(self, destHgDb, annotationType):
        self._needOptions("version")
        return "{}.{}.transMap{}.bigPsl".format(destHgDb, annotationType, self.version)

    def getMappedBigPslFileForDestHgDb(self, destHgDb, annotationType):
        return os.path.join(self.getMappedResultsForDestHgDb(destHgDb),
                            self._bigPslFileBasename(destHgDb, annotationType))

    def getMappedGbdbDir(self, destHgDb):
        self._needOptions("version")
        return os.path.join(self.gbdbDir, destHgDb, "transMap", self.version)

    def getMappedGbdbBigPslFile(self, destHgDb, annotationType):
        return os.path.join(self.getMappedGbdbDir(destHgDb),
                            self._bigPslFileBasename(destHgDb, annotationType))

    @property
    def mappedBigPslFile(self):
        self._needOptions("destHgDb", "annotationType")
        return self.getMappedBigPslFileForDestHgDb(self.destHgDb, self.annotationType)

    @property
    def mappingChainsDir(self):
        self._needOptions("srcHgDb", "destHgDb")
        return os.path.join(self.dataRootDir, "data/chains", self.srcHgDb, self.destHgDb)

    @property
    def mappingChains(self):
        self._needOptions("srcHgDb", "destHgDb", "chainType")
        chainFile = "{}.{}.{}.chain".format(self.srcHgDb, self.destHgDb, self.chainType)
        return os.path.join(self.mappingChainsDir, chainFile)

    @staticmethod
    def _twoBitSearch(twoBitPathPat, db):
        paths = [pat.format(db=db) for pat in twoBitPathPat]
        for twoBit in paths:
            if os.path.exists(twoBit):
                return twoBit
        raise Exception("can't find 2bit file for {db} in {paths}".format(db=db, paths=paths))

    @staticmethod
    def _chromSizesGet(chromSizesPat, db):
        path = chromSizesPat.format(db=db)
        if os.path.exists(path):
            return path
        raise Exception("can't find chrom size file for {db}: {path}".format(db=db, path=path))

    def getSrcTwoBit(self):
        "search for the src twobit file"
        self._needOptions("srcHgDb")
        return self._twoBitSearch(self.srcTwoBitPathPat, self.srcHgDb)

    def getSrcChromSizes(self):
        "get the src chrmo size file"
        self._needOptions("srcHgDb")
        return self._chromSizesGet(self.srcChromSizesPat, self.srcHgDb)

    def getDestTwoBit(self):
        "search for the dest twobit file"
        self._needOptions("destHgDb")
        return self._twoBitSearch(self.destTwoBitPathPat, self.destHgDb)

    def getDestChromSizes(self):
        "get the dest chrmo size file"
        self._needOptions("destHgDb")
        return self._chromSizesGet(self.destChromSizesPat, self.destHgDb)

    @property
    def mappingChainsIndexDb(self):
        return self.mappingChains() + ".db"

    @property
    def batchDir(self):
        self._needOptions("buildDir", "destHgDb", "srcHgDb", "annotationType")
        return os.path.join(self.buildDir, "batches", self.destHgDb, self.srcHgDb, str(self.annotationType))

    @property
    def batchParaDir(self):
        self._needOptions("batchGen")
        return os.path.join(self.batchDir, "para{}".format(self.batchGen))

    @property
    def batchParaFile(self):
        return os.path.join(self.batchParaDir, "transmap.batch")

    @property
    def batchResultsDir(self):
        return os.path.join(self.batchDir, "results")

    @property
    def numSeqsPerJob(self):
        "get the number of sequence per batch mapping job for this conf"
        if self.fixedNumSeqsPerJob is not None:
            return self.fixedNumSeqsPerJob
        self._needOptions("srcHgDb", "annotationType")
        if self.annotationType in (AnnotationType.rna, AnnotationType.est):
            return 2500
        else:
            return 500

    @property
    def jobPreBigPslDir(self):
        return os.path.join(self.batchResultsDir, "parts")

    def getJobPreBigPsl(self, startOid, endOid):
        return os.path.join(self.jobPreBigPslDir, "{}.{}.preBigPsl".format(startOid, endOid))

    def getBatchSrcHgDbPreBigPsl(self, srcHgDb):
        self._needOptions("destHgDb", "annotationType")
        return os.path.join(self.buildDir, "results", self.destHgDb,
                            "{}.{}.{}.preBigPsl".format(self.destHgDb, srcHgDb, self.annotationType))

    @property
    def batchPreBigPsl(self):
        self._needOptions("srcHgDb")
        return self.getBatchSrcHgDbPreBigPsl(self.srcHgDb)

    @property
    def bigTransMapAsPath(self):
        return os.path.expanduser("~/kent/src/hg/lib/bigTransMap.as")

    @property
    def statsDir(self):
        self._needOptions("dataRootDir")
        return os.path.normpath(os.path.join(self.dataRootDir, "stats"))

    @property
    def detailedStatsDir(self):
        return os.path.join(self.statsDir, "detailed")

    def getDetailedStatsTsv(self, destHgDb):
        return os.path.join(self.detailedStatsDir, "{}.detailed.tsv".format(destHgDb))

    @property
    def detailedStatsTsv(self):
        self._needOptions("destHgDb")
        return self.getDetailedStatsTsv(self.destHgDb)

    def _loadRequiredPreviousDestHgDbs(self):
        if self.prevDataRootDir is None:
            return frozenset()
        else:
            return frozenset(os.listdir(os.path.join(self.prevDataRootDir, "results/mapped")))

    @property
    def requiredPreviousDestHgDbs(self):
        if self.requiredPreviousDestHgDbsCache is None:
            self.requiredPreviousDestHgDbsCache = self._loadRequiredPreviousDestHgDbs()
        return self.requiredPreviousDestHgDbsCache


def transMapConfLoad(configPyFile, dataRootDir=None, srcHgDb=None, destHgDb=None,
                     annotationType=None, chainType=None):
    """
    Loads the configuration file which must have define a function getConfig,
    which takes the same signature as TransMapConf.__init__, with the exception
    of the version and batchGen, and then adds these parameters and returns an instance
    of the TransMapConf class."""
    getFuncKwargs = {"dataRootDir": dataRootDir,
                     "srcHgDb": srcHgDb,
                     "destHgDb": destHgDb,
                     "annotationType": annotationType,
                     "chainType": chainType}
    return evalConfigFunc(configPyFile, getFuncKwargs=getFuncKwargs)
