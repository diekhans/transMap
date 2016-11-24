"""
Database with range index into temporary mapping chains.  Chains are
not loaded into
"""
from collections import namedtuple
from pycbio.hgdata.hgLite import HgLiteTable
from pycbio.hgdata.rangeFinder import Binner


class MappingChainsDbTables(object):
    """tables names in a mappings sqlite3 database"""
    mappingChainsTbl = "mappingChains"


class MappingChainLoc(namedtuple("MappingChainsIndex",
                                 ("bin", "qName", "qStart", "qEnd", "offset", "length"))):
    """Index of chain record"""
    __slots__ = ()

    @staticmethod
    def factory(qName, qStart, qEnd, offset, length):
        return MappingChainLoc(Binner.calcBin(qStart, qEnd), qName, qStart, qEnd, offset, length)


class MappingChainsIndexDbTable(HgLiteTable):
    """
    Database table containing index into mapping chains.
    """
    __createSql = """CREATE TABLE {table} (
            bin int unsigned not null,
            qName text not null,
            qStart int unsigned not null,
            qEnd int unsigned not null,
            offset int unsigned not null,
            length int unsigned not null);"""
    __insertSql = """INSERT INTO {table} (bin, qName, qStart, qEnd, offset, length) VALUES (?, ?, ?, ?, ?, ?);"""
    __indexSql = """CREATE INDEX {table}_chrom_bin on {table} (qName, bin)"""

    def __init__(self, conn, table, create=False):
        super(MappingChainsIndexDbTable, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows, which should already have bin"""
        self._inserts(self.__insertSql, rows)

    def getRangeOverlap(self, qName, qStart, qEnd):
        """Get index entries for overlapping range"""
        sql = "SELECT * FROM {{table}} WHERE {};".format(Binner.getOverlappingSqlExpr("bin", "qName", "qStart", "qEnd", qName, qStart, qEnd))
        return self.queryRows(sql, lambda cur, row: MappingChainLoc(*row))
