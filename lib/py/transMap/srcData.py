"""
sqllite3 databases objects for transmap intermediate data.
"""

from collections import namedtuple
from transMap import alignIdToSrcId, srcIdToAccv, accvToAcc
from pycbio.hgdata.hgLite import HgLiteTable, SequenceDbTable, PslDbTable


# FIXME: maybe rename from *Lite to *Tables


class SourceDbTables(object):
    """tables names in a source data sqlite3 database"""
    srcAlignTbl = "srcAlign"
    srcXRefTbl = "srcXRef"
    srcMetadataTbl = "srcMetadata"
    srcSeqTbl = "srcSeq"


class TransMapSrcGene(namedtuple("TransMapSrcGene", ("srcId", "accv", "cds", "geneId", "geneName", "geneType", "transcriptType"))):
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
            geneId text,
            geneName text,
            geneType text,
            transcriptType text);"""
    __insertSql = """INSERT INTO {table} (srcId, accv, cds, geneId, geneName, geneType, transcriptType) VALUES (?, ?, ?, ?, ?, ?, ?);"""
    __indexSql = """CREATE UNIQUE INDEX {table}_srcId on {table} (srcId);"""

    def __init__(self, conn, table, create=False):
        super(TransMapSrcGeneLite, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or TransMapSrcGene objects"""
        self._inserts(self.__insertSql, rows)

    @staticmethod
    def loadStep(transMapSrcDbConn, metadataReader):
        "function to load and index"
        transMapGeneTbl = TransMapSrcGeneLite(transMapSrcDbConn, SourceDbTables.srcMetadataTbl, True)
        transMapGeneTbl.loads(list(metadataReader))
        transMapGeneTbl.index()


class TransMapSrcXRef(namedtuple("TransMapSrcXRef", ("srcAlignId", "srcId", "accv"))):
    """link between transmap source ids"""
    __slots__ = ()

    @staticmethod
    def fromSrcAlnId(srcAlignId):
        "construct from a srcAlignId"
        srcId = alignIdToSrcId(srcAlignId)
        return TransMapSrcXRef(srcAlignId, srcId, srcIdToAccv(srcId))


class TransMapSrcXRefLite(HgLiteTable):
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

    def __init__(self, conn, table, create=False):
        super(TransMapSrcXRefLite, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or TransMapSrcGene objects"""
        self._inserts(self.__insertSql, rows)

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


class TransMapSrcAlignLite(PslDbTable):
    """source alignments, in PSL format"""
    def __init__(self, conn, table, create=False):
        super(TransMapSrcAlignLite, self).__init__(conn, table, create)

    def getAllAccv(self):
        """get set of accv for PSLs that were loaded; use to restrict set for testing"""
        sql = """SELECT qName FROM {table};"""
        return frozenset([srcIdToAccv(alignIdToSrcId(row[0]))
                          for row in self.query(sql)])

    def getAllAcc(self):
        """get set of acc (no version) for PSLs that were loaded; use to restrict set for testing"""
        return frozenset([accvToAcc(accv) for accv in self.getAllAccv()])


def srcAlignXRefLoad(transMapSrcDbConn, alignReader):
    "load function the alignments and xrefs from a psl with srcAlignId in qName"
    psls = list(alignReader)
    srcXRefs = [TransMapSrcXRef.fromSrcAlnId(psl[9]) for psl in psls]

    with transMapSrcDbConn:
        srcAlignTbl = TransMapSrcAlignLite(transMapSrcDbConn, SourceDbTables.srcAlignTbl, True)
        srcAlignTbl.loads(psls)
        srcAlignTbl.index()

        srcXRefTbl = TransMapSrcXRefLite(transMapSrcDbConn, SourceDbTables.srcXRefTbl, True)
        srcXRefTbl.loads(srcXRefs)
        srcXRefTbl.index()


def getAccvSubselectClause(field, accvSet):
    return """({} in ({}))""".format(field, ",".join(['"{}"'.format(accv) for accv in accvSet]))


def loadSeqFa(tmpSeqFa, transMapSrcDbConn):
    seqTbl = SequenceDbTable(transMapSrcDbConn, SourceDbTables.srcSeqTbl, True)
    seqTbl.loadFastaFile(tmpSeqFa)
    seqTbl.index()


def querySrcPsls(transMapSrcDbConn):
    srcAlignTbl = TransMapSrcAlignLite(transMapSrcDbConn, SourceDbTables.srcAlignTbl)
    return srcAlignTbl.query("SELECT {} FROM {{table}};".format(TransMapSrcAlignLite.columnsNamesSql))
