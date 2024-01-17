"""
Database with range index into temporary mapping chains.  Chains are
not loaded into
"""
from collections import namedtuple
from pycbio.hgdata.hgSqlite import HgSqliteTable
from pycbio.hgdata.rangeFinder import calcBin, getOverlappingSqlExpr


class MappingChainsDbTables(object):
    """tables names in a mappings sqlite3 database"""
    mappingChainsTbl = "mappingChains"


class MappingChainLoc(namedtuple("MappingChainsIndex",
                                 ("bin", "qName", "qStart", "qEnd", "chainId", "offset", "length"))):
    """Index of chain record; we include bin just to simplify loading,
    query coordinates are on positive strand."""
    __slots__ = ()

    @staticmethod
    def factory(qName, qStart, qEnd, chainId, offset, length):
        "query coordinates must be positive strand."
        bin = calcBin(qStart, qEnd)
        return MappingChainLoc(bin, qName, qStart, qEnd, chainId, offset, length)


class MappingChainsIndexSqliteTable(HgSqliteTable):
    """
    Database table containing index into mapping chains.
    query confidantes are kept in terms of positive strand.
    """
    _createSql = """CREATE TABLE {table} (
            bin int unsigned not null,
            qName text not null,
            qStart int unsigned not null,
            qEnd int unsigned not null,
            chainId int unsigned not null,
            offset int unsigned not null,
            length int unsigned not null);"""
    _insertSql = """INSERT INTO {table} ({columns}) VALUES ({values});"""
    _indexSql = """CREATE INDEX {table}_chrom_bin on {table} (qName, bin)"""

    columnNames = ("bin", "qName", "qStart", "qEnd", "chainId", "offset", "length")

    def __init__(self, conn, table, create=False):
        super(MappingChainsIndexSqliteTable, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        self._create(self._createSql)

    def index(self):
        """create index after loading"""
        self._index(self._indexSql)

    def loads(self, rows):
        """load rows, which should already have bin"""
        self._inserts(self._insertSql, self.columnNames, rows)

    def getRangeOverlap(self, qName, qStart, qEnd):
        """Get index entries for overlapping range"""
        binWhere = getOverlappingSqlExpr("bin", "qName", "qStart", "qEnd", qName, qStart, qEnd)
        sql = "SELECT {{columns}} FROM {{table}} WHERE {};".format(binWhere)
        return self.queryRows(sql, self.columnNames, lambda cur, row: MappingChainLoc(*row))
