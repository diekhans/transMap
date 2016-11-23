""""Database with genomes and chains."""
import sqlite3
import re
from collections import namedtuple
from pycbio.sys.symEnum import SymEnum
from pycbio.sys import fileOps
from pycbio.hgdata.hgLite import HgLiteTable


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
        return AnnotationTypeSet([AnnotationType(v) for v in sval.split(',')])


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
                        ("srcDb", "destDb", "chainType", "chainFile", "netFile", "dist"))):
    "Describes the chains used in a mapping to a destDb"
    # FIXME make srcDb/destDb and chainsFinder.py queryHgDb targetHgDb consistent
    __slots__ = ()

    @staticmethod
    def rowFactory(cur, row):
        return Chains(*row)


class ChainsDbTable(HgLiteTable):
    """Interface todata on chains stored in the database"""
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
        super(ChainsDbTable, self).__init__(conn, table)
        if create:
            self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table"""
        self._inserts(self.__insertSql, rows)

    def queryByDbs(self, srcDb, destDb):
        sql = "SELECT * FROM {table} WHERE (srcDb=?) AND (destDb=?)"
        return self.queryRows(sql, Chains.rowFactory, srcDb, destDb)

    def queryByDbsType(self, srcDb, destDb, chainType):
        "return Chain or None"
        sql = "SELECT * FROM {table} WHERE (srcDb=?) AND (destDb=?) AND (chainType=?)"
        rows = list(self.queryRows(sql, Chains.rowFactory, srcDb, destDb, chainType))
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
        self.dbNum = self.dbParse(hgDb)[1]

    def __str__(self):
        return "{}\t{}\t{}\t{}\t".format(self.hgDb, self.clade, self.commonName, self.scientificName, str(self.annotationTypeSet))

    dbParseRe = re.compile("^([a-zA-Z]+)([0-9]+)$")

    @staticmethod
    def dbParse(dbName):
        "parse a genomedb name into (orgDbName, dbNum)"
        m = GenomeAsm.dbParseRe.match(dbName)
        if m is None:
            raise Exception("can't parse database name: " + dbName)
        return (m.group(1), int(m.group(2)))

    def isFinished(self):
        "is this genome considered finished"
        return self.commonName in ("Human", "Mouse")

    def matches(self, clades=None, annotationType=None):
        "does this match the specified filter sets"
        if (clades is not None) and (self.clade not in clades):
            return False
        if (annotationType is not None) and (len(self.annotationType & annotationType) == 0):
            return False
        return True


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


class Genomes(object):
    "object containing all genome definitions, loaded from database"
    def __init__(self):
        self.genomeAsms = dict()          # by hgDb
        self.orgs = dict()         # commmonName to Organism (-> GenomeDb)
        self.clades = set()        # all clades

    def addGenomeAsm(self, genomeAsm):
        "add a genome db"
        self.genomeAsms[genomeAsm.hgDb] = genomeAsm
        org = self.obtainOrganism(genomeAsm.commonName)
        org.add(genomeAsm)
        self.clades.add(genomeAsm.clade)

    def obtainOrganism(self, commonName):
        org = self.orgs.get(commonName)
        if org is None:
            org = self.orgs[commonName] = Organism(commonName)
        return org

    def finish(self):
        """Finish up, making all links"""
        for org in self.orgs.itervalues():
            org.sort()

        # paranoia
        self.clades = frozenset(self.clades)

    def getDbByName(self, hgDb):
        "lookup db name or None"
        return self.genomeAsms.get(hgDb)

    def dump(self, fh):
        "print for debugging purposes"
        fileOps.prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.itervalues():
            fileOps.prRowv(fh, "", org.commonName)
            for db in org.dbs:
                fileOps.prRowv(fh, "", "", db)
