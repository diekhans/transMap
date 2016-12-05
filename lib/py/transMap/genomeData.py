""""Database with genomes and chains."""
import sqlite3
from collections import namedtuple
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
                        ("srcHgDb", "destHgDb", "chainType", "chainFile", "netFile"))):
    "Describes the chains used in a mapping to a destHgDb"
    # FIXME make srcDb/destHgDb and chainsFinder.py queryHgDb targetHgDb consistent
    __slots__ = ()

    @staticmethod
    def rowFactory(cur, row):
        return Chains(*row)


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
        return GenomeAsm(*row)


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


class Genomes(object):
    """Object containing genome definitions, loaded from database for clades
    of interested. Also used to query relevant changes for mapping """
    def __init__(self, conf):
        self.conf = conf
        self.genomeAsms = dict()   # by hgDb
        self.orgs = dict()         # commmonName to Organism (-> GenomeDb)

        self.genomeDbConn = sqlite3.connect(conf.genomeDb)
        self.__loadGenomes(conf)
        self.__finish()
        self.chainsTbl = ChainsDbTable(self.genomeDbConn, GenomesDbTables.chainsTbl)

    def __loadGenomes(self):
        genomeAsmsTbl = GenomeAsmsDbTable(self.genomeDbConn, GenomesDbTables.genomeAsmsTbl)
        for genomeAsm in genomeAsmsTbl.queryByClades(self.genomeDbConn, self.conf.clades):
            self.__addGenomeAsm(genomeAsm)

    def __addGenomeAsm(self, genomeAsm):
        "add a genome db"
        self.genomeAsms[genomeAsm.hgDb] = genomeAsm
        org = self.__obtainOrganism(genomeAsm.commonName)
        org.add(genomeAsm)

    def __obtainOrganism(self, commonName):
        org = self.orgs.get(commonName)
        if org is None:
            org = self.orgs[commonName] = Organism(commonName)
        return org

    def __finish(self):
        """Finish up, making all links"""
        for org in self.orgs.itervalues():
            org.sort()

    def getDbByName(self, hgDb):
        "lookup db name or None"
        return self.genomeAsms.get(hgDb)

    def getChains(self, srcHgDb, destHgDb):
        return list(self.chainType.queryByDbs(srcHgDb, destHgDb))

    def dump(self, fh):
        "print for debugging purposes"
        fileOps.prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.itervalues():
            fileOps.prRowv(fh, "", org.commonName)
            for db in org.dbs:
                fileOps.prRowv(fh, "", "", db)
