import sys, os
myDir = os.path.normpath(os.path.dirname(__file__))
sys.path.append(myDir + "/../../../extern/pycbio/lib")

from transMap.GenomeDefs import CDnaTypes

class TransMap(object):
    "object that has common data and functions used my transMap rules"
    def __init__(self, exrun, genomeDefs):
        self.exrun = exrun
        self.genomeDefs = genomeDefs
        self.buildDir = "build/test"
        self.arch = os.uname()[4]
        
    def getDataDir(self, srcDb):
        return self.buildDir + "/data/" + srcDb

    def getSrcPsl(self, srcDb, cdnaType):
        path=self.getDataDir() + str(cdnaType) + ".psl.bz2"
        return self.exrun.getFile(path)
        
    def getSrcAlnStats(self, srcDb, cdnaType):
        path=self.getDataDir() + str(cdnaType) + ".pslstats.bz2"
        return self.exrun.getFile(path)
        
    def getSrcFa(self, srcDb, cdnaType):
        path=self.getDataDir() + str(cdnaType) + ".fa.bz2"
        return self.exrun.getFile(path)
        
    def getSrcCds(self, srcDb, cdnaType):
        path=self.getDataDir() + str(cdnaType) + ".cds.bz2"
        return self.exrun.getFile(path)
        
        
        
