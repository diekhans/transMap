""""Objects used to define genomes and what mappings can be performed.
The GenomeDefs object can be constructed and pickled by the genomeDefsMk
program to speed up loading this data.
"""
import re
from pycbio.sys.Enumeration import Enumeration
from pycbio.sys.MultiDict import MultiDict

# type of chains to use in doing mapping
ChainSelect = Enumeration("ChainSelect", ("all", "syn"))

# chain directions (rev must be swapped)
ChainDir = Enumeration("ChainDir", ("fwd", "rev"))

# cDNA set types
CDnaTypes = Enumeration("CDnaTypes", ("refSeq", "mrna", "splicedEst"))

class Chains(object):
    "Describes a set of chains used in a mapping to a destDb"

    def __init__(self, srcDb, destDb, chainSelect, chainDir):
        """srcDb and destDb are or will become GenomeDb objects, they can be
        initialized as names"""
        self.srcDb = srcDb
        self.destDb = destDb
        self.chainSelect = ChainSelect(str(chainSelect))
        self.chainDir = ChainDir(str(chainDir))

class GenomeDb(object):
    "All information about a genome database"

    dbParseRe = re.compile("^([a-zA-Z]+)([0-9]+)$")

    def __init__(self, db, clade, cdnaTypes):
        self.db = db
        m = self.dbParseRe.match(self.db)
        self.dbOrg = m.group(1)
        self.dbNum = m.group(2)
        self.clade = clade
        self.cdnaTypes = set()
        for cdna in cdnaTypes:
            self.cdnaTypes.add(CDnaTypes(str(cdna)))
        self.cdnaTypes = frozenset(self.cdnaTypes)
        self.srcDbChains = []  # chains to this database

class GenomeDefs(object):
    "object containing all genome definitions"
    def __init__(self):
        self.dbs = dict()
        self.orgs = MultiDict()    # GenomeDb by org, sorted newest to oldest
        self.clades = MultiDict()  # GenomeDb by clade
        self.chains = []           # Chains objects
        self.dirty = True

    def addGenome(self, gdb):
        "add a GenomeDb object"
        assert(isinstance(gdb, GenomeDb))
        self.dbs[gdb.db] = gdb
        self.orgs.add(gdb.dbOrg, gdb)
        self.dirty = True
        
    def addChains(self, chains):
        "add a Chains object"
        self.chains.append(chains)
        self.dirty = True
        
    def __mapDbName(self, db):
        "change string db name to GenomeDb, if needed"
        if isinstance(db, str):
            return self.dbs[db]
        else:
            return db

    def __finishChains(self):
        for ch in self.chains:
            ch.srcDb = self.__mapDbName(ch.srcDb)
            ch.desDb = self.__mapDbName(ch.destDb)
            ch.srcDb.srcDbChains.append(ch)

    def __getLatestDestDbs(self):
        dbs = set()
        for org in self.orgs.iterkeys():
            gdb = self.orgs[org][0]
            if len(gdb.srcDbChains) > 0:
                dbs.add(gdb)
        return dbs

    def __getDestDbs(self, destDbs):
        dbs = set()
        for gdb in map(self.__mapDbName, destDbs):
            if len(gdb.srcDbChains) == 0:
                raise Exception("no chains mapping to ", gbd.db)
            dbs.add(gdb)
        return dbs


    def finish(self, destDbs=None, requiredDestDbs=None):
        """Finish up and determine the set of destination genomes.  If not
        None, destDbs is a list dbs that should have genes mapped to the, The
        default is the lastest version of all with sources. If requiredDestDbs
        is not None, it's a list to add to the default.  Can be called again
        if object has more chains or genomes added."""

        # sort so newest db is first in each org list
        for org in self.orgs.iterkeys():
            self.orgs[org].sort(lambda a,b: -cmp(a.dbNum,b.dbNum))

        self.__finishChains()

        # build set of destDbs
        if destDbs != None:
            self.destDbs = self.__getDestDbs(destDbs)
        else:
            self.destDbs = self.__getLatestDestDbs()
        if requiredDestDbs != None:
            self.destDbs |= self.__getDestDbs(requiredDestDbs)
        self.dirty = False


