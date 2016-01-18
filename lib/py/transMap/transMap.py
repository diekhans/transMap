import os, sys
from .genomeDefs import AnnotationSet
from pycbio.sys.pipeline import Pipeline
from pycbio.sys import fileOps
from pycbio.tsv import TSVReader
from . import getTmpExt, runCmds


class TransMap(object):
    """Object that has common data and functions used my transMap rules.
    If exrun is None, then it returns string paths rather than references"""
    def __init__(self, exrun, genomeDefs, buildDir, mappings, clusterDir=None):
        self.exrun = exrun
        self.genomeDefs = genomeDefs
        self.buildDir = buildDir
        self.clusterDir = clusterDir
        self.mappings = mappings
        self.arch = os.uname()[4]

    minQSize = 20
    minCover = 0.20
    maxAligns = 100

    def __getFile(self, path):
        if self.exrun != None:
            return self.exrun.getFile(path)
        else:
            return path

    ###
    # source data
    ##
    def getDataDir(self, srcDb):
        return self.buildDir + "/data/" + srcDb.name

    def getDataPre(self, srcDb, annotationSet):
        return self.getDataDir(srcDb) + "/" + str(annotationSet)

    def getSrcPsl(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".psl.bz2")
        
    def getSrcSeqId(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".seqid.bz2")
        
    def getSrcAlnStats(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".pslstats.bz2")
        
    def getSrcFa(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".fa.bz2")
        
    def getSrcPolyA(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".polya.bz2")
        
    def getSrcMeta(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".meta.bz2")
        
    def getSrcCds(self, srcDb, annotationSet):
        return self.__getFile(self.getDataPre(srcDb, annotationSet) + ".cds.bz2")
        
    ###
    # cluster source data (partitioned by prefix).
    ##
    def getClusterDataDir(self, srcDb, annotationSet):
        return self.clusterDir + "/data/" + srcDb.name + "/" + str(annotationSet)

    def getClusterDataPre(self, srcDb, annotationSet, cdnaPart):
        return self.getClusterDataDir(srcDb, annotationSet) + "/" + cdnaPart

    def getClusterSrcPsl(self, srcDb, annotationSet, cdnaPart):
        return self.__getFile(self.getClusterDataPre(srcDb, annotationSet, cdnaPart) + ".psl")
        
    def getClusterSrcFa(self, srcDb, annotationSet, cdnaPart):
        return self.__getFile(self.getClusterDataPre(srcDb, annotationSet, cdnaPart) + ".fa")

    ##
    # intermediate transMap alignments before filtering
    ##
    def getClusterMappedDir(self, destDb, srcDb):
        return self.clusterDir + "/aligns/" + destDb.name + "/" + srcDb.name

    def getClusterMappedPre(self, destDb, srcDb, annotationSet, cdnaPart, genomePart=None):
        p = self.getClusterMappedDir(destDb, srcDb) + "/" + str(annotationSet) + "/" + cdnaPart
        if genomePart != None: 
            p += "/" + genomePart
        return p

    def getClusterMappedPsl(self, destDb, srcDb, annotationSet, cdnaPart, genomePart=None):
        return self.__getFile(self.getClusterMappedPre(destDb, srcDb, annotationSet, cdnaPart, genomePart) + ".psl")

    def getClusterMappedInfo(self, destDb, srcDb, annotationSet, cdnaPart, genomePart=None):
        return self.__getFile(self.getClusterMappedPre(destDb, srcDb, annotationSet, cdnaPart, genomePart) + ".mapinfo")

    ##
    # transMap alignments
    ##
    def getMappedDir(self, destDb, srcDb):
        return self.buildDir + "/aligns/" + destDb.name + "/" + srcDb.name

    def getMappedPre(self, destDb, srcDb, annotationSet):
        return self.getMappedDir(destDb, srcDb) + "/" + str(annotationSet)

    def getMappedPsl(self, destDb, srcDb, annotationSet):
        return self.__getFile(self.getMappedPre(destDb, srcDb, annotationSet) + ".psl.bz2")

    def getMappedInfo(self, destDb, srcDb, annotationSet):
        return self.__getFile(self.getMappedPre(destDb, srcDb, annotationSet) + ".mapinfo.bz2")


    ##
    # filtering options
    ##
    def getGlobalNearBest(self, destDb):
        "get globalNearBest cutoff"
        # finished targets are more strigent
        if destDb.isFinished():
            return 0.005
        else:
            return 0.01

    def getGlobalNearBestPreFilter(self, destDb):
        "get globalNearBest prefilter cutoff"
        # larger than final filter, since match stats columns not set
        return 2*self.getGlobalNearBest(destDb)

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
    "cache of Chroms, by GenomeDb"
    def __init__(self, dataDir):
        self.dataDir = dataDir

    def obtain(self, db):
        chroms = self.get(db)
        if chroms == None:
            chroms = self[db] = Chroms(self.dataDir + "/" + db.name + "/" + db.name + ".2bit")
        return chroms

class CDnaPartitions(object):
    def __init__(self, db, annotationSet, transMap, cdnaTargetSize):
        self.db = db
        self.annotationSet = annotationSet
        self.transMap = transMap
        self.cdnaTargetSize = cdnaTargetSize
        self.partList = transMap.getClusterDataDir(db, annotationSet) + "/parts.lst"
        if not os.path.exists(self.partList):
            self.__mkPartitions()
        self.parts = fileOps.readFileLines(self.partList)

    def __mkPartitions(self):
        outDir = self.transMap.getClusterDataDir(self.db, self.annotationSet)
        fileOps.ensureDir(outDir)
        partListTmp = self.partList + getTmpExt()
        runCmds(["pslFaPartition", "-partList="+partListTmp,
                 self.cdnaTargetSize,
                 self.transMap.getSrcPsl(self.db, self.annotationSet),
                 self.transMap.getSrcFa(self.db, self.annotationSet),
                 outDir])
        os.rename(partListTmp, self.partList)

class CDnaPartitionsCache(dict):
    "cache of CDnaPartitions, by (GenomeDb, AnnotationSet)"
    def __init__(self, transMap, cdnaTargetSize):
        self.transMap = transMap
        self.cdnaTargetSize = cdnaTargetSize

    def obtain(self, db, annotationSet):
        key = (db, annotationSet)
        cdnaParts = self.get(key)
        if cdnaParts == None:
            cdnaParts = self[key] = CDnaPartitions(db, annotationSet, self.transMap, self.cdnaTargetSize)
        return cdnaParts

    
