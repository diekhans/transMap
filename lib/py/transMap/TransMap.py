import os
from transMap.GenomeDefs import CDnaTypes

class TransMap(object):
    "object that has common data and functions used my transMap rules"
    def __init__(self, exrun, genomeDefs, buildDir, mappings):
        self.exrun = exrun
        self.genomeDefs = genomeDefs
        self.buildDir = buildDir
        self.mappings = mappings
        self.arch = os.uname()[4]
        
    def getDataDir(self, srcDb):
        return self.buildDir + "/data/" + srcDb.db

    def getDataPre(self, srcDb, cdnaType):
        return self.getDataDir(srcDb) + "/" + str(cdnaType)

    def getSrcPsl(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".psl.bz2"
        return self.exrun.getFile(path)
        
    def getSrcSeqId(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".seqid.bz2"
        return self.exrun.getFile(path)
        
    def getSrcAlnStats(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".pslstats.bz2"
        return self.exrun.getFile(path)
        
    def getSrcFa(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".fa.bz2"
        return self.exrun.getFile(path)
        
    def getSrcPolyA(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".polya.bz2"
        return self.exrun.getFile(path)
        
    def getSrcMeta(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".meta.bz2"
        return self.exrun.getFile(path)
        
    def getSrcCds(self, srcDb, cdnaType):
        path = self.getDataPre(srcDb, cdnaType) + ".cds.bz2"
        return self.exrun.getFile(path)
        
    def getMappedDir(self, destDb, srcDb):
        return self.buildDir + "/mapped/" + destDb.db + "/" + srcDb.db

    def getMappedPre(self, destDb, srcDb, cdnaType):
        return self.getMappedDir(destDb, srcDb) + "/" + str(cdnaType)

    def getMappedRawPsl(self, destDb, srcDb, cdnaType):
        path = self.getMappedPre(destDb, srcDb, cdnaType) + ".rawPsl.bz2"
        return self.exrun.getFile(path)

    def getMappedRawMapinfo(self, destDb, srcDb, cdnaType):
        path = self.getMappedPre(destDb, srcDb, cdnaType) + ".rawMapinfo.bz2"
        return self.exrun.getFile(path)

    def getMappedPsl(self, destDb, srcDb, cdnaType):
        path = self.getMappedPre(destDb, srcDb, cdnaType) + ".psl.bz2"
        return self.exrun.getFile(path)

    def getMappedFiltStats(self, destDb, srcDb, cdnaType):
        path = self.getMappedPre(destDb, srcDb, cdnaType) + ".filtStats"
        return self.exrun.getFile(path)

    def getMappedPslStats(self, destDb, srcDb, cdnaType):
        path = self.getMappedPre(destDb, srcDb, cdnaType) + ".pslstats.bz2"
        return self.exrun.getFile(path)

