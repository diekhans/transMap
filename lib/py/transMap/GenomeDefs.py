""""Objects used to define genomes and what mappings can be performed.
The GenomeDefs object can be constructed and pickled by the genomeDefsMk
program to speed up loading this data.
"""
import re
from pycbio.sys.Enumeration import Enumeration
from pycbio.sys.MultiDict import MultiDict
from pycbio.sys import setOps
from pycbio.sys.fileOps import prRowv


# type of chains to use in doing mapping
ChainSelect = Enumeration("ChainSelect", ("all", "syn"))

# chain directions (rev must be swapped)
ChainDir = Enumeration("ChainDir", ("fwd", "rev"))

# cDNA set types
CDnaTypes = Enumeration("CDnaTypes", ("refSeq", "mrna", "splicedEst"))

class Chains(object):
    "Describes a set of chains used in a mapping to a destDb"
    __slots__ = ("srcDb", "destDb", "dir")
    def __init__(self, srcDb, destDb, dir):
        """srcDb and destDb are or will become GenomeDb objects, they can be
        initialized as db names"""
        self.srcDb = srcDb
        self.destDb = destDb
        self.dir = ChainDir(str(dir))

    def __str__(self):
        if self.dir == ChainDir.fwd:
            return self.srcDb.db + " ==> " + self.destDb.db
        else:
            return self.srcDb.db + " <== " + self.destDb.db

    def __getstate__(self):
        return (self.srcDb, self.destDb, self.dir)

    def __setstate__(self, st):
        (self.srcDb, self.destDb, self.dir) = st

class GenomeDb(object):
    "All information about a genome database"

    dbParseRe = re.compile("^([a-zA-Z]+)([0-9]+)$")

    def __init__(self, db, clade, cdnaTypes):
        self.db = db
        m = self.dbParseRe.match(self.db)
        self.dbOrg = m.group(1)
        self.dbNum = m.group(2)
        self.clade = clade
        self.cdnaTypes = frozenset(cdnaTypes)
        self.srcDbChains = []  # chains to this database

    def __str__(self):
        return self.db + " " + self.clade + " " + setOps.setJoin(self.cdnaTypes,",")

class GenomeDefs(object):
    "object containing all genome definitions"
    def __init__(self):
        self.dbs = dict()
        self.orgs = MultiDict()    # GenomeDb by org, sorted newest to oldest
        self.chains = []           # all Chains objects
        self.srcDbs = set()        # all source dbs
        self.destDbs = set()       # all destination dbs

    def addGenomeDb(self, db, clade, cdnaTypes):
        "add a genome db"
        gdb = GenomeDb(db, clade,cdnaTypes)
        self.dbs[gdb.db] = gdb
        self.orgs.add(gdb.dbOrg, gdb)
        
    def addChains(self, srcDb, destDb, dir):
        "add a set of chains"
        self.chains.append(Chains(srcDb, destDb, dir))
        
    def __mapDbName(self, db):
        "change string db name to GenomeDb, if needed"
        if isinstance(db, str):
            return self.dbs[db]
        else:
            return db

    def __finishChains(self):
        for ch in self.chains:
            ch.srcDb = self.__mapDbName(ch.srcDb)
            ch.destDb = self.__mapDbName(ch.destDb)
            ch.srcDb.srcDbChains.append(ch)

    def finish(self):
        """Finish up, making all likes"""

        # sort so newest db is first in each org list
        for org in self.orgs.iterkeys():
            self.orgs[org].sort(lambda a,b: -cmp(a.dbNum,b.dbNum))

        self.__finishChains()
        for ch in self.chains:
            self.srcDbs.add(ch.srcDb)
            self.destDbs.add(ch.destDb)

    def __dbSetToStr(self, dbSet):
        "convert a set of GenomeDb to a string of db names"
        l = []
        for db in dbSet:
            l.append(db.db)
        l.sort()
        return ",".join(l)

    def dump(self, fh):
        "print for debugging purposes"
        prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.iterkeys():
            for db in self.orgs[org]:
                prRowv(fh, "", db)
        prRowv(fh, "GenomeDefs: chains:")
        for ch in self.chains:
            prRowv(fh, "", ch)
        prRowv(fh, "GenomeDefs: srcDbs: ", self.__dbSetToStr(self.srcDbs))
        prRowv(fh, "GenomeDefs: destDbs:", self.__dbSetToStr(self.destDbs))
        
