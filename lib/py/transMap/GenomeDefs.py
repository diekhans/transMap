""""Objects used to define genomes and what mappings can be performed.
The GenomeDefs object can be constructed and pickled by the genomeDefsMk
program to speed up loading this data.
"""
import os, re, cPickle, glob
from pycbio.sys.Enumeration import Enumeration
from pycbio.sys.MultiDict import MultiDict
from pycbio.sys import setOps,typeOps
from pycbio.sys.fileOps import prRowv

# FIXME: problem with Enumeration as variables rather than
# classes is serialization doesn't match!

# type of chains, which can either be existing or ones filtered on the fly.
ChainType = Enumeration("ChainType", ("all", "syn", "rbest"))

# cDNA set types
CDnaType = Enumeration("CDnaType", ("refSeq", "mrna", "splicedEst", "ucscGenes"))

def enumSetStr(eset):
    return "["+setOps.setJoin(eset, ",")+"]"

class Chains(object):
    "Describes a set of chains used in a mapping to a destDb"
    __slots__ = ("srcDb", "destDb", "chtype", "chainFile", "netFile", "dist")
    def __init__(self, srcDb, destDb, chtype, chainFile, netFile):
        """srcDb and destDb are GenomeDb objects"""
        assert(isinstance(srcDb, GenomeDb))
        assert(isinstance(destDb, GenomeDb))
        self.srcDb = srcDb
        self.destDb = destDb
        self.chtype = chtype
        self.chainFile = chainFile
        self.netFile = netFile
        self.dist = None

    def __str__(self):
        return self.srcDb.name + " ==> " + self.destDb.name + " [" + str(self.chtype) + "] |" + str(self.dist) + "|"

    def dump(self, fh):
        prRowv(fh, str(self), self.bzdir)

    def __getstate__(self):
        return (self.srcDb, self.destDb, self.chtype, self.chainFile, self.netFile, self.dist)

    def __setstate__(self, st):
        (self.srcDb, self.destDb, self.chtype, self.chainFile, self.netFile, self.dist) = st

class ChainsSet(object):
    """set of chains between two dbs.  There can be multiple chains
    due to different chain types."""
    
    # /cluster/data/hg18/bed/blastz.calJac1/axtChain/
    #     hg18.calJac1.all.chain.gz
    #     hg18.calJac1.net.gz
    #     hg18.calJac1.rbest.chain.gz
    #     hg18.calJac1.rbest.net.gz
    # /cluster/data/hg18/bed/blastz.ponAbe2/axtChain/
    #     hg18.ponAbe2.syn.net.gz
    
    __slots__ = ("srcDb", "destDb", "byType", "dist")
    def __init__(self, srcDb, destDb):
        self.srcDb = srcDb
        self.destDb = destDb
        self.byType = dict() # by types
        self.dist = None
        bzdir = self.__findChainsDir(self.srcDb, self.destDb)
        if bzdir != None:
            self.__addIfExists(bzdir, ChainType.all, "all.chain", "net")
            self.__addIfExists(bzdir, ChainType.syn, "all.chain", "syn.net")
            self.__addIfExists(bzdir, ChainType.rbest, "rbest.chain", "rbest.net")

    def __getstate__(self):
        return (self.srcDb, self.destDb, self.byType, self.dist)

    def __setstate__(self, st):
        (self.srcDb, self.destDb, self.byType, self.dist) = st

    def setDistance(self, dist):
        "set distance"
        self.dist = dist
        for ch in self.byType.itervalues():
            ch.dist = dist

    def __str__(self):
        return self.srcDb.name + " ==> " + self.destDb.name + " " + enumSetStr(set(self.byType.iterkeys())) + " |" + str(self.dist) + "|"

    @staticmethod
    def __findChainsDir(srcDb, destDb):
        "return path to chain directory if one can be found, otherwise None"
        base = "/cluster/data/" + destDb.name + "/bed/blastz." + srcDb.name
        if os.path.exists(base):
            return base
        p = base + ".swap"
        if os.path.exists(p):
            return p
        # look for data extensions w/wo .swap
        dirs = glob.glob(base + ".*-*-*")
        dirs.sort()
        if len(dirs) > 0:
            return dirs[-1]

        # check for missing symlink to names like blastzPanTro2.2006-03-28
        pat = "/cluster/data/" + destDb.name + "/bed/blastz" + srcDb.name[0].upper() + srcDb.name[1:] + ".*-*-*"
        dirs = glob.glob(pat)
        dirs.sort()
        if len(dirs) > 0:
            return dirs[-1]
        return None

    def __getBlastzFile(self, bzdir, ext):
        """generate and check for one of the blastz alignment files, with and
        without compressed extensions.  Return path if one exists, otherwise
        None"""
        bzpre = bzdir + "/axtChain/" + self.destDb.name + "." + self.srcDb.name + "."
        p = bzpre + ext + ".gz"  # common case
        if os.path.exists(p):
            return p
        p = bzpre + ext
        if os.path.exists(p):
            return p
        return None

    def __addIfExists(self, bzdir, chtype, chainExt, netExt):
        """add a set of nets/chains if they exist for the specified type"""
        chainFile = self.__getBlastzFile(bzdir, chainExt)
        if chainFile != None:
            netFile = self.__getBlastzFile(bzdir, netExt)
            if netFile != None:
                self.byType[chtype] = Chains(self.srcDb, self.destDb, chtype, chainFile, netFile)

    def getSupportingType(self, chtype):
        """get the chains that support the give type, or None.  This can return
        the all chains when then require more filtering"""
        chains = self.byType.get(chtype)
        if (chains == None) and (chtype == ChainType.syn):
            chains = self.byType[ChainType.all]
        return chains

    def haveChains(self):
        return len(self.byType) > 0

    def getDb(self, useSrcDb):
        "get db by useSrcDb"
        return self.srcDb if useSrcDb else self.destDb

    @staticmethod
    def distCmp(cs1, cs2):
        if (cs1.dist == None) and (cs2.dist == None):
            d = cmp(cs1.srcDb.name, cs2.srcDb.name)
            if d == 0:
                d = cmp(cs1.destDb.name, cs2.destDb.name)
            return d
        elif cs1.dist == None:
            return 1
        elif cs2.dist == None:
            return -1
        else:
            return cmp(cs1.dist, cs2.dist)

class OrgChainsSetMap(MultiDict):
    """map of Organism to list of ChainsSet, sorted newest to oldest, can be
    keyed and sorted by srcDb or destDb"""

    def __init__(self, db, srcDbIndexed):
        self.db = db
        self.srcDbIndexed = srcDbIndexed

    def __str__(self):
        if self.srcDbIndexed:
            l = [srcDb.name for chset in self]
        else:
            l = [destDb.name for chset in self]
        l.sort()
        return ",".join(l)

    def add(self, chainsSet):
        "add a new database"
        assert(chainsSet.getDb(not self.srcDbIndexed) == self.db)
        org = chainsSet.srcDb.org if self.srcDbIndexed else chainsSet.destDb.org
        MultiDict.add(self, org, chainsSet)

    def addAll(self, chainsSets):
        "add all chainsSets"
        for chset in chainsSets:
            self.add(chset)

    def sort(self):
        "sort all entries"
        for ent in self.iterentries():
            if self.srcDbIndexed:
                ent.sort(lambda a,b: -cmp(a.srcDb.dbNum,b.srcDb.dbNum))
            else:
                ent.sort(lambda a,b: -cmp(a.destDb.dbNum,b.destDb.dbNum))

    def getChainsSet(self, db):
        entries = self[db.org]
        if self.srcDbIndexed:
            for chset in entries:
                if chset.srcDb == db:
                    return chset
        else:
            for chset in entries:
                if chset.destDb == db:
                    return chset
        raise Exception("no chain from for " + str(db))

class GenomeDb(object):
    "All information about a genome database"

    dbParseRe = re.compile("^([a-zA-Z]+)([0-9]+)$")

    @staticmethod
    def dbParse(dbName):
        "parse a genomedb name into (orgDbName, dbNum)"
        m = GenomeDb.dbParseRe.match(dbName)
        if m == None:
            raise Exception("can't parse database name: " + dbName)
        return (m.group(1), int(m.group(2)))

    def __init__(self, dbName, clade, common, scientific, cdnaTypes=None):
        self.name = dbName
        self.orgDbName, self.dbNum = self.dbParse(self.name)
        self.clade = clade
        self.common = common
        self.scientific = scientific
        self.cdnaTypes = frozenset(cdnaTypes) if (cdnaTypes != None) else None
        # chains to/from this databases, list by orgDbName, sorted newest to oldest
        self.srcChainsSets = OrgChainsSetMap(self, True)
        self.destChainsSets = OrgChainsSetMap(self, False)

    def finish(self):
        "complete object"
        self.srcChainsSets.sort()
        self.destChainsSets.sort()

    def __str__(self):
        s = self.name + "\t" + self.clade + "\t'" + self.common + "'\t'" + self.scientific + "'\t"
        if self.cdnaTypes != None:
            s += setOps.setJoin(self.cdnaTypes,",")
        else:
            s += "[]"
        return s

    def isFinished(self):
        "is this genome considered finished"
        return (self.orgDbName == "hg") or  (self.orgDbName == "mm")

    def matches(self, clades=None, cdnaTypes=None):
        "does this match the specified filter sets"
        if (clades != None) and (self.clade not in clades):
            return False
        if (cdnaTypes != None) and (len(self.cdnaTypes & cdnaTypes) == 0):
            return False
        return True

class Organism(object):
    "databases for a given organism"

    def __init__(self, orgDbName):
        self.orgDbName = orgDbName
        self.dbs = []

    def add(self, db):
        "add a new database"
        assert(db.orgDbName == self.orgDbName)
        self.dbs.append(db)
        db.org = self

    def sort(self):
        "sort all entries, by reverse dbNum"
        self.dbs.sort(lambda a,b: -cmp(a.dbNum,b.dbNum))

class GenomeDefs(object):
    "object containing all genome definitions"
    def __init__(self):
        self.dbs = dict()
        self.orgs = dict()         # orgDbName to Organism (-> GenomeDb)
        self.chainsSets = []       # all ChainsSet objects
        self.srcChainsSets = set()        # all source dbs
        self.destChainsSets = set()       # all destination dbs
        self.clades = set()        # all clades
        self.cdnaTypes = set()     # all CDnaTypes

    def addGenomeDb(self, db):
        "add a genome db"
        self.dbs[db.name] = db
        self.obtainOrganism(db.orgDbName).add(db)
        self.clades.add(db.clade)
        self.cdnaTypes |= db.cdnaTypes

    def obtainOrganism(self, orgDbName):
        org = self.orgs.get(orgDbName)
        if org == None:
            org = self.orgs[orgDbName] = Organism(orgDbName)
        return org

    def buildChainSets(self):
        "build chain sets between all databases"
        for srcDb in self.dbs.itervalues():
            for destDb in self.dbs.itervalues():
                if srcDb.org != destDb.org:
                    chset = ChainsSet(srcDb, destDb)
                    if chset.haveChains():
                        srcDb.destChainsSets.add(chset)
                        destDb.srcChainsSets.add(chset)
                        self.chainsSets.append(chset)

    def finish(self):
        """Finish up, making all likes"""
        for org in self.orgs.itervalues():
            org.sort()
        for chset in self.chainsSets:
            self.srcChainsSets.add(chset.srcDb)
            self.destChainsSets.add(chset.destDb)
        for db in self.dbs.itervalues():
            db.finish()

        # paranoia
        self.srcChainsSets = frozenset(self.srcChainsSets)
        self.destChainsSets = frozenset(self.destChainsSets)
        self.clades = frozenset(self.clades)
        self.cdnaTypes = frozenset(self.cdnaTypes)

    def mapDbName(self, db):
        "change string db name to GenomeDb, if its not already a GenomeDb"
        if isinstance(db, str):
            return self.dbs[db]
        else:
            assert(isinstance(db, GenomeDb))
            return db

    def mapDbNames(self, dbs):
        """change list of db references, which can be either string names
        or GenomeDb objects, to a list of GenomeDb objects"""
        objs = []
        for db in dbs:
            objs.append(self.mapDbName(db))
        return objs

    def mapDbNamesSet(self, dbs):
        """change list of to db names or GenomeDb object, or None
        and return frozenset."""
        return frozenset(self.mapDbNames(typeOps.mkset(dbs)))

    def dump(self, fh):
        "print for debugging purposes"
        prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.iterkeys():
            for db in self.orgs[org].dbs:
                prRowv(fh, "", db, enumSetStr(db.cdnaTypes))
        prRowv(fh, "GenomeDefs: chains:")
        for chset in self.chainsSets:
            prRowv(fh, "", chset)

    def save(self, path):
        fh = open(path, "w")
        cPickle.dump((ChainType, CDnaType, self), fh)
        # FIXME , cPickle.HIGHEST_PROTOCOL causes error
        fh.close()

def load(path):
    fh = open(path)
    try:
        global ChainType, CDnaType
        ChainType, CDnaType, defs = cPickle.load(fh)
        return defs
    finally:
        fh.close()
