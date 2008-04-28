""""Objects used to define genomes and what mappings can be performed.
The GenomeDefs object can be constructed and pickled by the genomeDefsMk
program to speed up loading this data.
"""
import os, re, cPickle, glob
from pycbio.sys.Enumeration import Enumeration
from pycbio.sys.MultiDict import MultiDict
from pycbio.sys import setOps,typeOps
from pycbio.sys.fileOps import prRowv


# type of chains to use in doing mapping
ChainSelect = Enumeration("ChainSelect", ("all", "syn"))

# chain directions (rev must be swapped)
ChainOrient = Enumeration("ChainOrient", ("fwd", "rev"))

# cDNA set types
CDnaTypes = Enumeration("CDnaTypes", ("refSeq", "mrna", "splicedEst", "ucscGenes"))

class Chains(object):
    "Describes a set of chains used in a mapping to a destDb"
    __slots__ = ("srcDb", "destDb", "orient", "bzdir")
    def __init__(self, srcDb, destDb):
        """srcDb and destDb are GenomeDb objects"""
        assert(isinstance(srcDb, GenomeDb))
        assert(isinstance(destDb, GenomeDb))
        self.srcDb = srcDb
        self.destDb = destDb
        self.orient = None
        self.bzdir = None

    def __str__(self):
        if self.orient == None:
            sym = " ??? "
        elif self.orient == ChainOrient.fwd:
            sym = " ==> "
        else:
            sym = " <== "
        return self.srcDb.db + sym + self.destDb.db

    def dump(self, fh):
        prRowv(fh, str(self), self.bzdir)

    def __getstate__(self):
        return (self.srcDb, self.destDb, self.orient, self.bzdir)

    def __setstate__(self, st):
        (self.srcDb, self.destDb, self.orient, self.bzdir) = st

    def __findChainsDir(self, srcDb, destDb):
        "return path to chain directory if one can be found, otherwise None"
        base = "/cluster/data/" + destDb.db + "/bed/blastz." + srcDb.db
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

        # ditch, for missign symlink to names like blastzPanTro2.2006-03-28
        pat = "/cluster/data/" + destDb.db + "/bed/blastz" + srcDb.db[0].upper() + srcDb.db[1:] + ".*-*-*"
        dirs = glob.glob(pat)
        dirs.sort()
        if len(dirs) > 0:
            return dirs[-1]
        return None

    def setup(self):
        """determine the chain directory and orientation by looking for the directories.  Set
        the fields in object and return True or False"""
        self.bzdir = self.__findChainsDir(self.srcDb, self.destDb)
        if self.bzdir != None:
            self.orient = ChainOrient.fwd
        else:
            self.bzdir = self.__findChainsDir(self.destDb, self.srcDb)
            if self.bzdir != None:
                self.orient = ChainOrient.rev
        return (self.bzdir != None)
        
    def __getBlastzPre(self):
        if self.orient == ChainOrient.fwd:
            p = self.bzdir + "/axtChain/" + self.destDb.db + "." + self.srcDb.db
        else:
            p = self.bzdir + "/axtChain/" + self.srcDb.db + "." + self.destDb.db
        return p

    def getChainFile(self):
        f = self.__getBlastzPre() + ".all.chain.gz"
        if os.path.exists(f):
            return f
        else:
            return self.__getBlastzPre() + ".all.chain"

    def haveChainFile(self):
        return os.path.exists(self.getChainFile())
        
    def getNetFile(self):
        f = self.__getBlastzPre() + ".net.gz"
        if os.path.exists(f):
            return f
        else:
            return self.__getBlastzPre() + ".net"

    def haveNetFile(self):
        return os.path.exists(self.getNetFile())

    def haveChainNet(self):
        return self.haveChainFile() and self.haveNetFile()

    def sortCmp(a, b):
        d = cmp(a.srcDb.dbOrg, b.srcDb.dbOrg)
        if d == 0:
            d = cmp(a.srcDb.dbNum, b.srcDb.dbNum)
            if d == 0:
                d = cmp(a.destDb.dbOrg, b.destDb.dbOrg)
                if d == 0:
                    d = cmp(a.destDb.dbNum, b.destDb.dbNum)
        return d
        
        
    sortCmp = staticmethod(sortCmp)
        
class OrgDbMap(MultiDict):
    "map of dbOrg to list of GenomeDb, sorted newest to oldest"

    def add(self, db):
        "add a new database"
        MultiDict.add(self, db.dbOrg, db)

    def sort(self):
        "sort all entries"
        for ent in self.iterentries():
            ent.sort(lambda a,b: -cmp(a.dbNum,b.dbNum))

class OrgDbChainMap(MultiDict):
    "map of dbOrg to list of Chains, sorted newest to oldest"

    def __init__(self, srcDbIndexed):
        self.srcDbIndexed = srcDbIndexed

    def add(self, chains):
        "add a new database"
        db = chains.srcDb if self.srcDbIndexed else chains.destDb
        MultiDict.add(self, db.dbOrg, chains)

    def sort(self):
        "sort all entries"
        for ent in self.iterentries():
            if self.srcDbIndexed:
                ent.sort(lambda a,b: -cmp(a.srcDb.dbNum,b.srcDb.dbNum))
            else:
                ent.sort(lambda a,b: -cmp(a.destDb.dbNum,b.destDb.dbNum))

class GenomeDb(object):
    "All information about a genome database"

    dbParseRe = re.compile("^([a-zA-Z]+)([0-9]+)$")

    def __init__(self, db, clade, cdnaTypes=None):
        self.db = db
        m = self.dbParseRe.match(self.db)
        if m == None:
            raise Exception("can't parse database name: " + db)
        self.dbOrg = m.group(1)
        self.dbNum = int(m.group(2))
        self.clade = clade
        self.cdnaTypes = frozenset(cdnaTypes) if (cdnaTypes != None) else None
        # chains to/from this databases, list by dbOrg, sorted newest to oldest
        self.srcDbs = OrgDbChainMap(True)
        self.destDbs = OrgDbChainMap(False)

    def finish(self):
        "complete object"
        self.srcDbs.sort()
        self.destDbs.sort()

    def __str__(self):
        s = self.db + " " + self.clade
        if self.cdnaTypes != None:
            s + " " + setOps.setJoin(self.cdnaTypes,",")
        return s

    def isFinished(self):
        "is this genome considered finished"
        return (self.dbOrg == "hg") or  (self.dbOrg == "mm")

    def matches(self, clades=None, cdnaTypes=None):
        "does this match the specified filter sets"
        if (clades != None) and (self.clade not in clades):
            return False
        if (cdnaTypes != None) and (len(self.cdnaTypes & cdnaTypes) == 0):
            return False
        return True

class GenomeMappings(list):
    "list of chains, along with sets of databases, to map"
    def __init__(self):
        self.srcDbs = set()
        self.destDbs = set()
        self.chainsMap = set()  # unordered, but quick lookup, self is sorted

    def _addMapping(self, chains):
        "add a mapping"
        self.append(chains)
        self.chainsMap.add(chains)
        self.srcDbs.add(chains.srcDb)
        self.destDbs.add(chains.destDb)

    def finish(self):
        "sort by srcDb, destDb"
        self.sort(Chains.sortCmp)

    def dump(self, fh):
        "print for debugging purposes"
        prRowv(fh, "Mappings:")
        for ch in self:
            fh.write("\t")
            fh.write(str(ch))
            if not ch.haveChainFile():
                fh.write("\tmissingChain: " + ch.getChainFile())
            if not ch.haveNetFile():
                fh.write("\tmissingNet: " + ch.getNetFile())
            fh.write("\n")

class GenomeDefs(object):
    "object containing all genome definitions"
    def __init__(self):
        self.dbs = dict()
        self.orgs = OrgDbMap()     # GenomeDb by org, sorted newest to oldest
        self.chains = []           # all Chains objects
        self.srcDbs = set()        # all source dbs
        self.destDbs = set()       # all destination dbs
        self.clades = set()        # all clades
        self.cdnaTypes = set()     # all CDnaTypes 

    def addGenomeDb(self, db):
        "add a genome db"
        self.dbs[db.db] = db
        self.orgs.add(db)
        self.clades.add(db.clade)
        self.cdnaTypes |= db.cdnaTypes
        
    def addChains(self, chains):
        "add a set of chains"
        self.chains.append(chains)
        
    def __mapDbName(self, db):
        "change string db name to GenomeDb, if needed"
        if isinstance(db, str):
            return self.dbs[db]
        else:
            return db

    def __mapDbNames(self, dbs):
        """change list of db references, which can be either string names
        or GenomeDb objects, to a list of GenomeDb objects"""
        objs = []
        for db in dbs:
            objs.append(self.__mapDbName(db))
        return objs

    def __mapDbNamesSet(self, dbs):
        """change list of to db names or GenomeDb object, or None
        and return frozenset."""
        return frozenset(self.__mapDbNames(typeOps.mkset(dbs)))

    def __finishChains(self):
        for ch in self.chains:
            ch.srcDb = self.__mapDbName(ch.srcDb)
            ch.destDb = self.__mapDbName(ch.destDb)
            ch.srcDb.destDbs.add(ch)
            ch.destDb.srcDbs.add(ch)
        self.chains.sort(Chains.sortCmp)

    def finish(self):
        """Finish up, making all likes"""

        self.orgs.sort()
        self.__finishChains()
        for ch in self.chains:
            self.srcDbs.add(ch.srcDb)
            self.destDbs.add(ch.destDb)
        for db in self.dbs.itervalues():
            db.finish()

        # paranoia
        self.srcDbs = frozenset(self.srcDbs)
        self.destDbs = frozenset(self.destDbs)
        self.clades = frozenset(self.clades)
        self.cdnaTypes = frozenset(self.cdnaTypes)

    def __dbSetToStr(self, dbSet):
        "convert a set of GenomeDb to a string of db names"
        l = []
        for db in dbSet:
            l.append(db.db)
        l.sort()
        return ",".join(l)

    def __useMapping(self, chains, clades, cdnaTypes, skipDbs, inclMissingChains, inclRevChains):
        return ((chains.srcDb not in skipDbs)
                and (inclRevChains or (chains.orient == ChainOrient.fwd))
                and (inclMissingChains or (chains.haveChainFile() and chains.haveNetFile()))
                and chains.srcDb.matches(clades, cdnaTypes))

    def __addDbMappings(self, gms, destDb, clades, cdnaTypes, skipDbs, inclMissingChains, inclRevChains):
        "add a GenomeMappings for a destDb"
        # add newest srcDb
        for ent in destDb.srcDbs.iterentries():
            if self.__useMapping(ent[0], clades, cdnaTypes, skipDbs, inclMissingChains, inclRevChains):
                gms._addMapping(ent[0])

    def __getDestDbs(self, onlyDestDbs, otherDestDbs, clades, skipDbs):
        "get set of destDbs, given restrictions"
        destDbs = set(otherDestDbs)
        if onlyDestDbs != None:
            destDbs |= onlyDestDbs
        else:
            for orgDbs in self.orgs.iterentries():
                if orgDbs[0].clade in clades:
                    destDbs.add(orgDbs[0])
        destDbs = destDbs.difference(skipDbs)
        return destDbs

    def getMappingSets(self, onlyDestDbs=None, otherDestDbs=None, clades=None, cdnaTypes=None,
                       skipDbs=None, inclMissingChains=False, inclRevChains=False):
        """Get a GenomeMappings object with chains describing possible
        mappings given the restrictions.  Only the latest db for an organism
        is used unless an earlier db is specified in otherDestDbs"""
        if onlyDestDbs != None:
            onlyDestDbs = self.__mapDbNamesSet(onlyDestDbs)
        # default to make code simpler
        otherDestDbs = self.__mapDbNamesSet(otherDestDbs)
        skipDbs = self.__mapDbNamesSet(skipDbs)
        clades = setOps.mkfzset(clades) if (clades != None) else self.clades
        cdnaTypes = setOps.mkfzset(cdnaTypes) if (cdnaTypes != None) else self.cdnaTypes

        gms = GenomeMappings()
        for destDb in self.__getDestDbs(onlyDestDbs, otherDestDbs, clades, skipDbs):
            self.__addDbMappings(gms, destDb, clades, cdnaTypes, skipDbs, inclMissingChains, inclRevChains)
        gms.finish()
        return gms

    def dump(self, fh):
        "print for debugging purposes"
        prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.iterkeys():
            for db in self.orgs[org]:
                prRowv(fh, "", db)
        prRowv(fh, "GenomeDefs: chains:")
        for ch in self.chains:
            prRowv(fh, "", ch, ch.bzdir)
        prRowv(fh, "GenomeDefs: srcDbs: ", self.__dbSetToStr(self.srcDbs))
        prRowv(fh, "GenomeDefs: destDbs:", self.__dbSetToStr(self.destDbs))

def load(path):
    fh = open(path)
    try:
        return cPickle.load(fh)
    finally:
        fh.close()
