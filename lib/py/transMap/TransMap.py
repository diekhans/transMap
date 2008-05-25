import os
from transMap.GenomeDefs import CDnaTypes
from pycbio.sys.Pipeline import Pipeline
from pycbio.tsv.TSVReader import TSVReader

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

class Chroms(list):
    """object contain list of chromosomes and sizes, sorted by descending size"""
    def __init__(self, genome2bit):
        pl = Pipeline([["twoBitInfo", genome2bit, "stdout"], ["fgrep", "-v", "_hap"], ["sort","-k 2,2nr"]])
        for row in TSVReader(genome2bit, columns=("chrom", "size"), typeMap={"size":int}, inFh=pl):
            self.append(row)

    def getRange(self, firstChrom, lastChrom):
        "get chroms in sort list bounds by firstChrom and lastChrom"
        l = len(self)
        i = 0
        while (i < l) and (self[i].chrom != firstChrom):
            i += 1
        if i >= l:
            raise Exception(firstChrom + " not found")
        chroms = []
        while i < l:
            chroms.append(self[i].chrom)
            if self[i].chrom == lastChrom:
                break
            i += 1
        if i >= l:
            raise Exception(lastChrom + " not found after " + firstChrom)
        return chroms

    def partition(self, targetSize, maxPerJob):
        """return list of (start, end); targetSize is used to combine scaffolds,
        maxPerJob prevents command lines from overflowing"""
        parts = []
        l = len(self)
        i = 0;
        while i < l:
            start = end = self[i].chrom
            size = self[i].size
            j = i+1
            cnt = 0
            while (j < l) and ((size+self[j].size) < targetSize) and (cnt < maxPerJob):
                end = self[j].chrom
                size += self[j].size
                j += 1
                cnt += 1
            parts.append((start, end))
            i = j
        return parts

class ChromsCache(dict):
    "cache of Chroms, by genome"
    def __init__(self, dataDir):
        self.dataDir = dataDir

    def obtain(self, db):
        chroms = self.get(db)
        if chroms == None:
            chroms = self[db] = Chroms(self.dataDir + "/" + db.db + "/" + db.db + ".2bit")
        return chroms
