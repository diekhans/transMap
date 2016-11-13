""""Objects used to define genomes and what mappings can be performed.
The GenomeDefs object can be constructed and pickled by the genomeDefsMk
program to speed up loading this data.
"""
import os
import sqlite3
import glob
from collections import namedtuple
from pycbio.sys.symEnum import SymEnum
from pycbio.sys.multiDict import MultiDict
from pycbio.sys import setOps, typeOps, fileOps
from pycbio.hgdata.hgLite import HgLiteTable


class GenomesDbTables(object):
    """tables names in a genomes sqlite3 database"""
    chainsTbl = "chains"


class ChainType(SymEnum):
    "type of chains, which can either be existing or ones filtered on the fly."
    all = 1
    syn = 2
    rbest = 3


sqlite3.register_adapter(ChainType, lambda chainType: str(chainType))
sqlite3.register_converter("chainType", ChainType)


class AnnSetType(SymEnum):
    "annotation set type"
    refseq = 1
    rna = 2
    est = 3
    gencode = 4
    ensembl = 5


class Chains(namedtuple("Chains",
                        ("srcDb", "destDb", "chainType", "chainFile", "netFile", "dist"))):
    "Describes the chains used in a mapping to a destDb"
    # FIXME make srcDb/destDb and chainsFinder.py queryHgDb targetHgDb consistent
    __slots__ = ()

    def __str__(self):
        return self.srcDb + " ==> " + self.destDb + " [" + str(self.chainType) + "] |" + str(self.dist) + "|"


class ChainsSet(object):
    """set of chains between two dbs.  There can be multiple chains
    due to different chain types."""

    # /hive/data/genomes/hg18/bed/blastz.calJac1/axtChain/
    #     hg18.calJac1.all.chain.gz
    #     hg18.calJac1.net.gz
    #     hg18.calJac1.rbest.chain.gz
    #     hg18.calJac1.rbest.net.gz
    # /hive/data/genomes/hg18/bed/blastz.ponAbe2/axtChain/
    #     hg18.ponAbe2.syn.net.gz

    __slots__ = ("srcDb", "destDb", "byType", "dist")

    def __init__(self, srcDb, destDb):
        self.srcDb = srcDb
        self.destDb = destDb
        self.byType = dict()  # by types
        self.dist = None
        blastzDir = self.__findChainsDir()
        if blastzDir is not None:
            self.__addIfExists(blastzDir, ChainType.all, "all.chain", "net")
            self.__addIfExists(blastzDir, ChainType.syn, "all.chain", "syn.net")
            self.__addIfExists(blastzDir, ChainType.rbest, "rbest.chain", "rbest.net")

    def __getstate__(self):
        return (self.srcDb, self.destDb, self.byType, self.dist)

    def __setstate__(self, st):
        (self.srcDb, self.destDb, self.byType, self.dist) = st

    def setDistance(self, dist):
        "set distance"
        self.dist = dist
        for ch in self.byType.itervalues():
            ch.dist = dist

    def getTypesStr(self):
        "string representation of chain types that are available"
        return setOps.setJoin(set(self.byType.iterkeys()))

    def __str__(self):
        return self.srcDb.name + " ==> " + self.destDb.name + " [" + self.getTypesStr() + "] |" + str(self.dist) + "|"

    def __chainsNetsExist(self, blastzDir):
        """ make sure both all.chain.gz and net.gz exist, as there was a case
        with only the rbest without an all.chain in one case"""
        return ((self.__getBlastzFile(blastzDir, "all.chain") is not None)
                and (self.__getBlastzFile(blastzDir, "net") is not None))

    def __findChainsDirMeth(self, meth):
        "Return path to chain directory for an alignment method."
        base = "/hive/data/genomes/" + self.destDb.name + "/bed/" + meth + "." + self.srcDb.name
        if os.path.exists(base):
            return base
        blastzDir = base + ".swap"
        if self.__chainsNetsExist(blastzDir):
            return blastzDir
        # look for data extensions w/wo .swap
        blastzDirs = glob.glob(base + ".*-*-*")
        blastzDirs.sort(reverse=True)
        for blastzDir in blastzDirs:
            if self.__chainsNetsExist(blastzDir):
                return blastzDir

        # check for missing symlink to names like blastzPanTro2.2006-03-28
        pat = "/hive/data/genomes/" + self.destDb.name + "/bed/" + meth + self.srcDb.name[0].upper() + self.srcDb.name[1:] + ".*-*-*"
        blastzDirs = glob.glob(pat)
        blastzDirs.sort()
        for blastzDir in blastzDirs:
            if self.__chainsNetsExist(blastzDir):
                return blastzDir
        return None

    def __findChainsDir(self):
        "return path to chain directory if one can be found, otherwise None"
        p = self.__findChainsDirMeth("lastz")
        if p is not None:
            p = self.__findChainsDirMeth("blastz")
        return p

    def __getBlastzFile(self, blastzDir, ext):
        """generate and check for one of the blastz alignment files, with and
        without compressed extensions.  Return path if one exists, otherwise
        None"""
        bzpre = blastzDir + "/axtChain/" + self.destDb.name + "." + self.srcDb.name + "."
        p = bzpre + ext + ".gz"  # common case
        if os.path.exists(p):
            return p
        p = bzpre + ext
        if os.path.exists(p):
            return p
        return None

    def __addIfExists(self, blastzDir, chainType, chainExt, netExt):
        """add a set of nets/chains if they exist for the specified type"""
        chainFile = self.__getBlastzFile(blastzDir, chainExt)
        if chainFile is not None:
            netFile = self.__getBlastzFile(blastzDir, netExt)
            if netFile is not None:
                self.byType[chainType] = Chains(self.srcDb, self.destDb, chainType, chainFile, netFile)

    def getSupportingType(self, chainType):
        """get the chains that support the give type, or None.  This can return
        the all chains when then require more filtering"""
        assert(isinstance(chainType, ChainType))
        chains = self.byType.get(chainType)
        if (chains is None) and (chainType == ChainType.syn):
            chains = self.byType[ChainType.all]
        return chains

    def haveChains(self):
        return len(self.byType) > 0

    def getDb(self, useSrcDb):
        "get db by useSrcDb"
        return self.srcDb if useSrcDb else self.destDb

    @staticmethod
    def distCmp(cs1, cs2):
        if (cs1.dist is None) and (cs2.dist is None):
            d = cmp(cs1.srcDb.name, cs2.srcDb.name)
            if d == 0:
                d = cmp(cs1.destDb.name, cs2.destDb.name)
            return d
        elif cs1.dist is None:
            return 1
        elif cs2.dist is None:
            return -1
        else:
            return cmp(cs1.dist, cs2.dist)

    def dump(self, fh, prefix=""):
        if len(self.byType) == 0:
            fh.write(prefix + self.srcDb.name + " ==> " + self.destDb.name + " [None]\n")
        else:
            for chains in self.byType.itervalues():
                fh.write(prefix + str(chains) + "\n")


class ChainsTbl(HgLiteTable):
    """Interface and cache of data on chains stored in the database"""
    __createSql = """CREATE TABLE {table} (
            srcDb text not null,
            destDb text not null,
            chainType text not null,
            chainFile text not null,
            netFile text not null,
            dist float);"""
    __insertSql = """INSERT INTO {table} (srcDb, destDb, chainType, chainFile, netFile, dist) VALUES (?, ?, ?, ?, ?, ?);"""
    __indexSql = ["""CREATE INDEX {table}_srcDb on {table} (srcDb);""",
                  """CREATE INDEX {table}_destDb on {table} (destDb);"""]

    def __init__(self, conn, table, create=False):
        super(ChainsTbl, self).__init__(conn, table)
        if create:
            self._create(self.__createSql)
        self.allChains = []

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table"""
        self._inserts(self.__insertSql, rows)


class ChainsSetMap(MultiDict):
    """map of Organism to list of ChainsSet, sorted newest to oldest, can be
    keyed and sorted by srcDb or destDb"""

    def __init__(self, db, srcDbIndexed):
        self.db = db
        self.srcDbIndexed = srcDbIndexed

    def __str__(self):
        if self.srcDbIndexed:
            l = [chset.srcDb.name for chset in self.itervalues()]
        else:
            l = [chset.destDb.name for chset in self.itervalues()]
        l.sort()
        return ",".join(l)

    def addChainSet(self, chainsSet):
        "add a new database"
        assert(chainsSet.getDb(not self.srcDbIndexed) == self.db)
        org = chainsSet.srcDb.org if self.srcDbIndexed else chainsSet.destDb.org
        self.add(org, chainsSet)

    def addChainSets(self, chainsSets):
        "add all chainsSets"
        for chset in chainsSets:
            self.addChainSet(chset)

    def sort(self):
        "sort all entries"
        for ent in self.iterentries():
            if self.srcDbIndexed:
                ent.sort(lambda a, b: -cmp(a.srcDb. dbNum, b.srcDb.dbNum))
            else:
                ent.sort(lambda a, b: -cmp(a.destDb. dbNum, b.destDb.dbNum))

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

    def dump(self, fh):
        fileOps.prRowv(fh, "ChainsSetMap", self.db)
        for chainsSet in self.itervalues():
            chainsSet.dump(fh, "\t")


class GenomeDb(object):
    "All information about a genome database"

    def __init__(self, dbName, org, clade, commonName, scientificName, annSetType=None):
        self.name = dbName
        self.org = org
        self.dbNum = self.dbParse(dbName)[1]
        self.clade = clade
        self.commonName = commonName
        self.scientificName = scientificName
        self.annSetType = frozenset(annSetType) if (annSetType is not None) else None
        # chains to/from this databases, list by Organism, sorted newest to oldest
        # self.srcChainsSets = OrgChainsSetMap(self, True)
        # self.destChainsSets = OrgChainsSetMap(self, False)

    def finish(self):
        "complete object"
        self.srcChainsSets.sort()
        self.destChainsSets.sort()

    def __str__(self):
        s = self.name + "\t" + self.clade + "\t'" + self.commonName + "'\t'" + self.scientificName + "'\t"
        if self.annSetType is not None:
            s += setOps.setJoin(self.annSetType, ",")
        else:
            s += "[]"
        return s

    def isFinished(self):
        "is this genome considered finished"
        return self.org.commonName in ("Human", "Mouse")

    def matches(self, clades=None, annSetType=None):
        "does this match the specified filter sets"
        if (clades is not None) and (self.clade not in clades):
            return False
        if (annSetType is not None) and (len(self.annSetType & annSetType) == 0):
            return False
        return True

    def __dumpChainsSets(self, fh, label, orgCsMap):
        fileOps.prLine(fh, "\t", label)
        for cs in orgCsMap.itervalues():
            fileOps.prLine(fh, "\t\t", cs)

    def dump(self, fh):
        fileOps.prLine(fh, str(self))
        self.__dumpChainsSets(fh, "srcChainsSets", self.srcChainsSets)
        self.__dumpChainsSets(fh, "destChainsSets", self.destChainsSets)


class Organism(object):
    """Databases for a given organism. We use common name to tie assembly dbs together, since db prefix
    has changed between releases: Baboon papHam1 ->papAnu2
    """

    def __init__(self, commonName):
        self.commonName = commonName
        self.dbs = []

    def add(self, db):
        "add a new database"
        assert(db.org == self)
        self.dbs.append(db)

    def sort(self):
        "sort all entries, by reverse dbNum"
        self.dbs.sort(lambda a, b: -cmp(a.dbNum, b.dbNum))


class GenomeDefs(object):
    "object containing all genome definitions"
    def __init__(self):
        self.genomes = dict()          # by name
        self.orgs = dict()         # commmonName to Organism (-> GenomeDb)
        self.clades = set()        # all clades
        self.annSetType = set()     # all AnnSetTypes

    def addGenomeDb(self, dbName, clade, commonName, scientificName, annSetType):
        "add a genome db"
        org = self.obtainOrganism(commonName)
        db = GenomeDb(dbName, org, clade, commonName, scientificName, annSetType)
        self.dbs[db.name] = db
        org.add(db)
        self.clades.add(db.clade)
        self.annSetType |= db.annSetType

    def obtainOrganism(self, commonName):
        org = self.orgs.get(commonName)
        if org is None:
            org = self.orgs[commonName] = Organism(commonName)
        return org

    def buildChainSets(self):
        "build chain sets between all databases"
        for srcDb in self.dbs.itervalues():
            for destDb in self.dbs.itervalues():
                if srcDb.org != destDb.org:
                    chset = ChainsSet(srcDb, destDb)
                    if chset.haveChains():
                        srcDb.destChainsSets.addChainSet(chset)
                        destDb.srcChainsSets.addChainSet(chset)

    def finish(self):
        """Finish up, making all likes"""
        for org in self.orgs.itervalues():
            org.sort()
        for db in self.dbs.itervalues():
            db.finish()

        # paranoia
        self.clades = frozenset(self.clades)
        self.annSetType = frozenset(self.annSetType)

    def getDbByName(self, dbName):
        "lookup db name or None"
        return self.dbs.get(dbName)

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
        """take list of db names or GenomeDb object, or None
        and return frozenset."""
        return frozenset(self.mapDbNames(typeOps.mkset(dbs)))

    def dump(self, fh):
        "print for debugging purposes"
        fileOps.prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.itervalues():
            fileOps.prRowv(fh, "", org.commonName)
            for db in org.dbs:
                fileOps.prRowv(fh, "", "", db)
