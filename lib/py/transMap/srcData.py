"""
sqllite3 databases objects for transmap intermediate data.
"""

from collections import namedtuple
from transMap import alignIdToSrcId, srcIdToAccv, accvToAcc
from pycbio.hgdata.hgLite import HgLiteTable, SequenceDbTable, PslDbTable


# FIXME: just use standard table names here.


class SourceDbTables(object):
    """tables names in a source data sqlite3 database"""
    srcAlignTbl = "srcAlign"
    srcXRefTbl = "srcXRef"
    srcMetadataTbl = "srcMetadata"
    srcSeqTbl = "srcSeq"


class SrcMetadata(namedtuple("SrcMetadata", ("srcId", "accv", "cds", "geneName", "geneId", "geneType", "transcriptType"))):
    """metedata on gene or mRNA """
    __slots__ = ()


class SrcMetadataDbTable(HgLiteTable):
    """
    Source metadata database table
    """
    __createSql = """CREATE TABLE {table} (
            srcId text not null,
            accv text not null,
            cds text,
            geneName text,
            geneId text,
            geneType text,
            transcriptType text);"""
    __insertSql = """INSERT INTO {table} ({columns}) VALUES ({values});"""
    __indexSql = """CREATE UNIQUE INDEX {table}_srcId on {table} (srcId);"""

    columnNames = ("srcId", "accv", "cds", "geneName", "geneId", "geneType",
                   "transcriptType")

    def __init__(self, conn, table, create=False):
        super(SrcMetadataDbTable, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or SrcMetadata objects"""
        self._inserts(self.__insertSql, self.columnNames, rows)

    def getBySrcId(self, srcId):
        "return row or error if not found"
        sql = """SELECT {columns} FROM {table} WHERE srcId = ?"""
        rows = list(self.queryRows(sql, self.columnNames, lambda cur, row: SrcMetadata(*row), srcId))
        if len(rows) == 0:
            raise Exception("can't find metadata for {}".format(srcId))
        return rows[0]

    @staticmethod
    def loadStep(srcDbConn, metadataReader):
        "function to load and index"
        srcMetadataTbl = SrcMetadataDbTable(srcDbConn, SourceDbTables.srcMetadataTbl, True)
        srcMetadataTbl.loads(list(metadataReader))
        srcMetadataTbl.index()


class SrcXRef(namedtuple("SrcXRef", ("srcAlignId", "srcId", "accv"))):
    """link between transmap source ids"""
    __slots__ = ()

    @staticmethod
    def fromSrcAlignId(srcAlignId):
        "construct from a srcAlignId"
        srcId = alignIdToSrcId(srcAlignId)
        return SrcXRef(srcAlignId, srcId, srcIdToAccv(srcId))


class SrcXRefDbTable(HgLiteTable):
    """
    transmap source id and accession
    """
    __createSql = """CREATE TABLE {table} (
            srcAlignId text not null,
            srcId text not null,
            accv text);"""
    __insertSql = """INSERT INTO {table} (srcAlignId, srcId, accv) VALUES (?, ?, ?);"""
    __indexSql = ["""CREATE UNIQUE INDEX {table}_srcAlignId on {table} (srcAlignId);""",
                  """CREATE INDEX {table}_accv on {table} (accv);"""]
    columnNames = ("srcAlignId", "srcId", "accv")

    def __init__(self, conn, table, create=False):
        super(SrcXRefDbTable, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or SrcMetadata objects"""
        self._inserts(self.__insertSql, self.columnNames, rows)

    def getSrcIds(self):
        "get generator over unique source ids"
        sql = "SELECT distinct(srcId) FROM {table};"
        for row in self.query(sql, []):
            yield row[0]

    def getAccvs(self):
        "get generator over unique source accv"
        sql = "SELECT distinct(accv) FROM {table};"
        for row in self.query(sql, []):
            yield row[0]


class SrcAlignDbTable(PslDbTable):
    """Source alignments, in PSL format.  Normally sorted by target
    to speed up alignment pipeline"""
    def __init__(self, conn, table, create=False):
        super(SrcAlignDbTable, self).__init__(conn, table, create)

    def getAllAccv(self):
        """get set of accv for PSLs that were loaded; use to restrict set for testing"""
        sql = """SELECT qName FROM {table};"""
        return frozenset([srcIdToAccv(alignIdToSrcId(row[0]))
                          for row in self.query(sql, [])])

    def getAllAcc(self):
        """get set of acc (no version) for PSLs that were loaded; use to restrict set for testing"""
        return frozenset([accvToAcc(accv) for accv in self.getAllAccv()])


def srcAlignXRefLoad(srcDbConn, alignReader):
    "load function the alignments and xrefs from a psl with srcAlignId in qName"
    psls = list(alignReader)
    srcXRefs = [SrcXRef.fromSrcAlignId(psl[9]) for psl in psls]

    with srcDbConn:
        srcAlignTbl = SrcAlignDbTable(srcDbConn, SourceDbTables.srcAlignTbl, True)
        srcAlignTbl.loads(psls)
        srcAlignTbl.index()

        srcXRefTbl = SrcXRefDbTable(srcDbConn, SourceDbTables.srcXRefTbl, True)
        srcXRefTbl.loads(srcXRefs)
        srcXRefTbl.index()


def getAccvSubselectClause(field, accvSet):
    return """({} in ({}))""".format(field, ",".join(['"{}"'.format(accv) for accv in accvSet]))


def loadSeqFa(tmpSeqFa, srcDbConn):
    seqTbl = SequenceDbTable(srcDbConn, SourceDbTables.srcSeqTbl, True)
    seqTbl.loadFastaFile(tmpSeqFa)
    seqTbl.index()


def querySrcPsls(srcDbConn):
    srcAlignTbl = SrcAlignDbTable(srcDbConn, SourceDbTables.srcAlignTbl)
    return srcAlignTbl.query("SELECT {columns} FROM {table};", SrcAlignDbTable.columnNames)
