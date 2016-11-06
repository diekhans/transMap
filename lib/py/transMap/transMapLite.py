"""
sqllite3 databases objects for transmap intermediate data.
"""

from collections import namedtuple
from pycbio.hgdata.hgLite import HgLiteTable, PslLite
from transMap import parseTransMapId


class SourceDbTables(object):
    """tables names in a source dataset sqlite3 database"""
    srcAlignsTbl = "srcAligns"
    srcMetaDataTbl = "srcMetaData"
    srcSeqTbl = "srcSeqs"


class TransMapGene(namedtuple("TransMapGene", ("id", "cds", "db", "geneName", "geneId"))):
    """metedata on gene or mRNA """
    __slots__ = ()


class TransMapGeneLite(HgLiteTable):
    """
    Source gene information
    """
    __createSql = """CREATE TABLE {table} (
            id text not null,
            cds text,
            db text not null,
            geneName text,
            geneId text);"""
    __insertSql = """INSERT INTO {table} (id, cds, db, geneName, geneId) VALUES (?, ?, ?, ?, ?);"""
    __indexSql = """CREATE UNIQUE INDEX {table}_id on {table} (id);"""

    def __init__(self, conn, table, create=False):
        super(TransMapGeneLite, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or TransMapGene objects"""
        self.inserts(self.__insertSql, rows)

    def index(self):
        """create index after loading"""
        self.execute(self.__indexSql)


def getTransMapSrcAccs(srcDbConn):
    """get the transmap accession for all source accessions"""
    alignSrcTbl = PslLite(srcDbConn, SourceDbTables.srcAlignsTbl)
    return frozenset(alignSrcTbl.queryRows("SELECT qName FROM {table}",
                                           lambda cur, row: parseTransMapId(row[0])[1]))
