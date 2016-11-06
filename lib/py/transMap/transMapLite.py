"""
sqllite3 databases objects for transmap intermediate data.
"""

from collections import namedtuple
from pycbio.hgdata.hgLite import HgLiteTable


class SourceDbTables(object):
    """tables names in a source dataset sqlite3 database"""
    srcAlignsTbl = "srcAligns"
    srcXRefTbl = "srcXRefs"
    srcMetaDataTbl = "srcMetaData"
    srcSeqTbl = "srcSeqs"


class TransMapSrcGene(namedtuple("TransMapSrcGene", ("srcId", "accv", "cds", "geneName", "geneId"))):
    """metedata on gene or mRNA """
    __slots__ = ()


class TransMapSrcGeneLite(HgLiteTable):
    """
    Source gene information
    """
    __createSql = """CREATE TABLE {table} (
            srcId text not null,
            accv text not null,
            cds text,
            geneName text,
            geneId text);"""
    __insertSql = """INSERT INTO {table} (srcId, accv, cds, geneName, geneId) VALUES (?, ?, ?, ?, ?);"""
    __indexSql = """CREATE UNIQUE INDEX {table}_srcId on {table} (srcId);"""

    def __init__(self, conn, table, create=False):
        super(TransMapSrcGeneLite, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or TransMapSrcGene objects"""
        self.inserts(self.__insertSql, rows)

    def index(self):
        """create index after loading"""
        self.execute(self.__indexSql)


class TransMapSrcXRef(namedtuple("TransMapSrcXRef", ("srcAlnId", "srcId", "accv"))):
    """link between transmap source ids"""
    __slots__ = ()


class TransMapSrcXRefLite(HgLiteTable):
    """
    transmap source id and accession
    """
    __createSql = """CREATE TABLE {table} (
            srcAlnId text not null,
            srcId text not null,
            accv text);"""
    __insertSql = """INSERT INTO {table} (srcAlnId, srcId, accv) VALUES (?, ?, ?);"""
    __indexSql = ["""CREATE UNIQUE INDEX {table}_srcAlnId on {table} (srcAlnId);""",
                  """CREATE INDEX {table}_accv on {table} (accv);"""]

    def __init__(self, conn, table, create=False):
        super(TransMapSrcXRefLite, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or TransMapSrcGene objects"""
        self.inserts(self.__insertSql, rows)

    def index(self):
        """create index after loading"""
        self.executes(self.__indexSql)

    def getSrcIds(self):
        "get generator over unique source ids"
        sql = "SELECT distinct(srcId) FROM {table};"
        for row in self.query(sql):
            yield row[0]

    def getAccvs(self):
        "get generator over unique source accv"
        sql = "SELECT distinct(accv) FROM {table};"
        for row in self.query(sql):
            yield row[0]
