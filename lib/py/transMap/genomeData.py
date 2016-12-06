""""Database with genomes and chains."""
import sqlite3
from collections import namedtuple, defaultdict
from pycbio.sys.symEnum import SymEnum
from pycbio.sys import fileOps
from pycbio.hgdata.hgLite import HgLiteTable
from transMap import hgDbNameParse


class GenomesDbTables(object):
    """tables names in a genomes sqlite3 database"""
    chainsTbl = "chains"
    genomeAsmsTbl = "genomeAsms"


class AnnotationType(SymEnum):
    "annotation type for a genome assembly"
    refseq = 0x01
    rna = 0x02
    est = 0x04
    gencode = 0x08
    ensembl = 0x10


class AnnotationTypeSet(frozenset):
    "annotation types set for a genome assembly"

    def __str__(self):
        values = []
        for value in self:
            values.append(str(value))
        return ",".join(values)

    @staticmethod
    def fromStr(sval):
        if sval is None:
            return AnnotationTypeSet()
        else:
            return AnnotationTypeSet([AnnotationType(str(v)) for v in sval.split(',')])


# FIXME register_converter doesn't seem to be working due to using rowFactory.  Pick one approach
sqlite3.register_adapter(AnnotationTypeSet, str)
sqlite3.register_converter("annotationTypes", AnnotationTypeSet.fromStr)


class ChainType(SymEnum):
    "type of chains, which can either be existing or ones filtered on the fly."
    all = 1
    syn = 2
    rbest = 3


sqlite3.register_adapter(ChainType, lambda chainType: str(chainType))
sqlite3.register_converter("chainType", ChainType)


class Chains(namedtuple("Chains",
                        ("srcHgDb", "destHgDb", "chainType", "chainFile", "netFile"))):
    "Describes the chains used in a mapping to a destHgDb"
    # FIXME make srcDb/destHgDb and chainsFinder.py queryHgDb targetHgDb consistent
    __slots__ = ()

    @staticmethod
    def rowFactory(cur, row):
        return Chains(row[0], row[1], ChainType(str(row[2])), row[3], row[4])

    def compatibleChainType(self, chainType):
        # syn is made from all
        return (((chainType == ChainType.syn) and (self.chainType == ChainType.all))
                or (chainType == self.chainType))

class ChainsDbTable(HgLiteTable):
    """Interface todata on chains stored in the database"""
    __createSql = """CREATE TABLE {table} (
            srcHgDb text not null,
            destHgDb text not null,
            chainType text not null,
            chainFile text not null,
            netFile text not null);"""
    __insertSql = """INSERT INTO {table} (srcHgDb, destHgDb, chainType, chainFile, netFile) VALUES (?, ?, ?, ?, ?);"""
    __indexSql = ["""CREATE INDEX {table}_srcHgDb on {table} (srcHgDb);""",
                  """CREATE INDEX {table}_destHgDb on {table} (destHgDb);"""]

    def __init__(self, conn, table, create=False):
        super(ChainsDbTable, self).__init__(conn, table)
        if create:
            self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table"""
        self._inserts(self.__insertSql, rows)

    def queryByDbs(self, srcHgDb, destHgDb):
        sql = "SELECT * FROM {table} WHERE (srcHgDb=?) AND (destHgDb=?)"
        return self.queryRows(sql, Chains.rowFactory, srcHgDb, destHgDb)

    def queryByDestHgDb(self, destHgDb):
        sql = "SELECT * FROM {table} WHERE (destHgDb=?)"
        return self.queryRows(sql, Chains.rowFactory, destHgDb)

    def queryByDbsType(self, srcHgDb, destHgDb, chainType):
        "return Chain or None"
        sql = "SELECT * FROM {table} WHERE (srcHgDb=?) AND (destHgDb=?) AND (chainType=?)"
        rows = list(self.queryRows(sql, Chains.rowFactory, srcHgDb, destHgDb, chainType))
        if len(rows) == 0:
            return None
        elif len(rows) == 1:
            return rows[0]
        else:
            raise Exception("too many rows returned")


class GenomeAsm(namedtuple("GenomeDb",
                           ("hgDb", "clade", "commonName", "scientificName", "annotationTypeSet"))):
    "Information about a genome assembly, stored in genome database "

    def __init__(self, hgDb, clade, commonName, scientificName, annotationTypeSet):
        super(GenomeAsm, self).__init__(hgDb, clade, commonName, scientificName, annotationTypeSet)
        self.dbNum = hgDbNameParse(hgDb)[1]

    def __str__(self):
        return "{}\t{}\t{}\t{}\t".format(self.hgDb, self.clade, self.commonName, self.scientificName, str(self.annotationTypeSet))

    def isFinished(self):
        "is this genome considered finished"
        return self.commonName in ("Human", "Mouse")

    @staticmethod
    def rowFactory(cur, row):
        return GenomeAsm(row[0], row[1], row[2], row[3], AnnotationTypeSet.fromStr(row[4]))


class GenomeAsmsDbTable(HgLiteTable):
    """Interface todata on chains stored in the database"""
    __createSql = """CREATE TABLE {table} (
            hgDb text not null,
            clade text not null,
            commonName text not null,
            scientificName text not null,
            annotationTypeSet text);"""
    __insertSql = """INSERT INTO {table} (hgDb, clade, commonName, scientificName, annotationTypeSet) VALUES (?, ?, ?, ?, ?);"""
    __indexSql = ["""CREATE INDEX {table}_hgDb on {table} (hgDb);""",
                  """CREATE INDEX {table}_commonName on {table} (commonName);""",
                  """CREATE INDEX {table}_annotationTypeSet on {table} (annotationTypeSet);"""]

    def __init__(self, conn, table, create=False):
        super(GenomeAsmsDbTable, self).__init__(conn, table)
        if create:
            self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table"""
        self._inserts(self.__insertSql, rows)

    def queryByClades(self, clades):
        sql = "SELECT * FROM {{table}} WHERE clade in ({})".format(",".join((len(clades) * ["?"])))
        return self.queryRows(sql, GenomeAsm.rowFactory, *clades)


class Organism(object):
    """Databases for a given organism. We use common name to tie assembly dbs together, since db prefix
    has changed between releases: Baboon papHam1 ->papAnu2
    """

    def __init__(self, commonName):
        self.commonName = commonName
        self.genomeAsms = []

    def add(self, genomeAsm):
        "add a genome assembly"
        self.genomeAsms.append(genomeAsm)

    def sort(self):
        "sort all entries, by reverse dbNum"
        self.genomeAsms.sort(lambda a, b: -cmp(a.dbNum, b.dbNum))


class Mapping(namedtuple("Mapping",
                         ("destHgDb", "srcHgDb", "annotationType", "chains"))):
    """Mapping to a destHgDb from a srcHgDb for an annotation type"""
    pass


class Genomes(object):
    """Object containing genome definitions, loaded from database for clades
    of interested. Also used to query relevant changes for mapping """
    def __init__(self, conf):
        self.conf = conf
        self.genomeAsms = dict()   # by hgDb
        self.orgs = dict()         # commmonName to Organism (-> GenomeDb)
        self.hgDbToOrg = dict()

        self.genomeDbConn = sqlite3.connect(conf.genomeDb)
        self.__loadGenomes()
        self.__finish()
        self.chainsTbl = ChainsDbTable(self.genomeDbConn, GenomesDbTables.chainsTbl)

    def __loadGenomes(self):
        genomeAsmsTbl = GenomeAsmsDbTable(self.genomeDbConn, GenomesDbTables.genomeAsmsTbl)
        for genomeAsm in genomeAsmsTbl.queryByClades(self.conf.clades):
            self.__addGenomeAsm(genomeAsm)

    def __addGenomeAsm(self, genomeAsm):
        "add a genome db"
        self.genomeAsms[genomeAsm.hgDb] = genomeAsm
        org = self.__obtainOrganism(genomeAsm.commonName)
        org.add(genomeAsm)
        self.hgDbToOrg[genomeAsm.hgDb] = org

    def __obtainOrganism(self, commonName):
        org = self.orgs.get(commonName)
        if org is None:
            org = self.orgs[commonName] = Organism(commonName)
        return org

    def __finish(self):
        """Finish up, making all links"""
        for org in self.orgs.itervalues():
            org.sort()

    def getAsmByName(self, hgDb):
        "lookup db name or None"
        return self.genomeAsms.get(hgDb)

    def findCurrentDestHgDbs(self):
        """find list of current destHgDb.  This is the newest from the target clades,
        plus previous ones required in the config.
        """
        currentDestHgDbs = set()
        for org in self.orgs.itervalues():
            currentDestHgDbs.add(org.genomeAsms[0].hgDb)
        for hgDb in self.conf.requiredPreviousDestHgDbs:
            if self.getAsmByName(hgDb) is None:
                raise Exception("config requiredPreviousDestHgDbs {} is not a valid hgDb".format(hgDb))
            currentDestHgDbs.add(hgDb)
        return currentDestHgDbs

    def __getDestDbChains(self, destHgDb):
        return self.chainsTbl.queryByDestHgDb(destHgDb)

    def __getChainsBySrcOrg(self, destHgDb):
        # skip srcHgDb not in table, these are outside of the specified clades
        chainsBySrcOrg = defaultdict(list)
        for chain in self.__getDestDbChains(destHgDb):
            if chain.srcHgDb in self.hgDbToOrg:
                chainsBySrcOrg[self.hgDbToOrg[chain.srcHgDb]].append(chain)
        return chainsBySrcOrg

    def __getPreferedChainTypes(self, destHgDb, srcHgDb):
        destClade = self.genomeAsms[destHgDb].clade
        srcClade = self.genomeAsms[srcHgDb].clade
        return self.conf.cladePreferedChains[frozenset([destClade, srcClade])]

    def __getDestHgDbSrcHgPreferedChain(self, destHgDb, srcHgDb, srcChainsList):
        # get most appropriate chain for the clade differences
        print destHgDb, srcHgDb,  self.__getPreferedChainTypes(destHgDb, srcHgDb)
        print "\n".join([str(c) for c in srcChainsList])
        for preferedChainType in self.__getPreferedChainTypes(destHgDb, srcHgDb):
            hit = [srcChain for srcChain in srcChainsList if srcChain.compatibleChainType(preferedChainType)]
            if len(hit) > 0:
                return hit[0]
        raise Exception("can't find preferedChainType for {} ({}) to {} ({}) in {}"
                        .format(destHgDb, self.genomeAsms[destHgDb].clade,
                                srcHgDb, self.genomeAsms[srcHgDb].clade,
                                srcChainsList))

    def __getDestHgDbSrcHgMappings(self, destHgDb, srcHgDb, srcChainsList):
        preferedChain = self.__getDestHgDbSrcHgPreferedChain(destHgDb, srcHgDb, srcChainsList)
        # build for each annotation set in src
        srcAsm = self.getAsmByName(srcHgDb)
        for annotationType in srcAsm.annotationTypeSet:
            yield Mapping(destHgDb, srcHgDb, annotationType, preferedChain)

    def __getDestHgDbSrcOrgMappings(self, destHgDb, srcOrg, srcOrgChainsList):
        srcHgDb = srcOrg.genomeAsms[0].hgDb # newest
        srcHgDbChains = [c for c in srcOrgChainsList  if c.srcHgDb == srcHgDb]
        if len(srcHgDbChains) > 0:
            return list(self.__getDestHgDbSrcHgMappings(destHgDb, srcHgDb, srcHgDbChains))
        else:
            return []

    def __getDestHgDbMappings(self, destHgDb):
        destHgDbMappings = []
        chainsBySrcOrg = self.__getChainsBySrcOrg(destHgDb)
        for srcOrg in chainsBySrcOrg.iterkeys():
            destHgDbMappings.extend(self.__getDestHgDbSrcOrgMappings(destHgDb, srcOrg, chainsBySrcOrg[srcOrg]))
        return destHgDbMappings

    def getCurrentMappings(self):
        currentMappings = []
        for destHgDb in self.findCurrentDestHgDbs():
            currentMappings.extend(self.__getDestHgDbMappings(destHgDb))
        return currentMappings

    def dump(self, fh):
        "print for debugging purposes"
        fileOps.prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.itervalues():
            fileOps.prRowv(fh, "", org.commonName)
            for db in org.dbs:
                fileOps.prRowv(fh, "", "", db)
