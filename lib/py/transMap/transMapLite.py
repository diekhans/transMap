"""
sqllite3 databases objects for transmap intermediate data.
"""

from collections import namedtuple
from pycbio.hgdata.hgLite import HgLiteTable
from transMap import alignIdToSrcId, srcIdToAccv, accvToAcc
from pycbio.hgdata.hgLite import SequenceLite, PslLite
from pycbio.sys import fileOps


class SourceDbTables(object):
    """tables names in a source dataset sqlite3 database"""
    srcAlignsTbl = "srcAligns"
    srcXRefTbl = "srcXRefs"
    srcMetadataTbl = "srcMetadata"
    srcSeqTbl = "srcSeqs"


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
            transcriptType test);"""
    __insertSql = """INSERT INTO {table} (srcId, accv, cds, geneId, geneName, geneType, transcriptType) VALUES (?, ?, ?, ?, ?, ?, ?);"""
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

    @staticmethod
    def fromSrcAlnId(srcAlnId):
        "construct from a srcAlnId"
        srcId = alignIdToSrcId(srcAlnId)
        return TransMapSrcXRef(srcAlnId, srcId, srcIdToAccv(srcId))


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


def loadAlignsXRefs(transMapSrcDbConn, alignReader):
    "load the alignments and xrefs from a psl with srcAlnId in qName"
    psls = list(alignReader)
    srcXRefs = [TransMapSrcXRef.fromSrcAlnId(psl[9]) for psl in psls]

    srcAlignsTbl = PslLite(transMapSrcDbConn, SourceDbTables.srcAlignsTbl, create=True)
    srcAlignsTbl.loads(psls)
    srcAlignsTbl.index()

    srcXRefTbl = TransMapSrcXRefLite(transMapSrcDbConn, SourceDbTables.srcXRefTbl, True)
    srcXRefTbl.loads(srcXRefs)
    srcXRefTbl.index()


def loadSrcGeneMetadata(transMapSrcDbConn, metadataReader):
    transMapGeneTbl = TransMapSrcGeneLite(transMapSrcDbConn, SourceDbTables.srcMetadataTbl, True)
    transMapGeneTbl.loads(list(metadataReader))
    transMapGeneTbl.index()


def getSrcAlignAccv(transMapSrcDbConn):
    """get set of accv for PSLs that were loaded; use to restrict set for testing"""
    srcAlignsTbl = PslLite(transMapSrcDbConn, SourceDbTables.srcAlignsTbl)
    sql = """SELECT qName FROM {table};"""
    return frozenset([srcIdToAccv(alignIdToSrcId(row[0]))
                      for row in srcAlignsTbl.query(sql)])

def getSrcAlignAcc(transMapSrcDbConn):
    """get set of acc (no version) for PSLs that were loaded; use to restrict set for testing"""
    return frozenset([accvToAcc(accv) for accv in getSrcAlignAccv(transMapSrcDbConn)])

def getAccvSubselectClause(field, accvSet):
    return """({} in ({}))""".format(field, ",".join(['"{}"'.format(accv) for accv in accvSet]))


def getSrcAccs(transMapSrcDbConn, accvFh):
    for accv in TransMapSrcXRefLite(transMapSrcDbConn, SourceDbTables.srcXRefTbl).getAccvs():
        accvFh.write("{}\n".format(accv))


def getSrcAccsFile(annSetType, transMapSrcDbConn):
    tmpAccvFile = fileOps.tmpFileGet("transMap.{}.".format(annSetType), ".acc")
    with open(tmpAccvFile, "w") as accvFh:
        getSrcAccs(transMapSrcDbConn, accvFh)
    return tmpAccvFile


def loadSeqFa(tmpSeqFa, transMapSrcDbConn):
    seqTbl = SequenceLite(transMapSrcDbConn, SourceDbTables.srcSeqTbl, True)
    seqTbl.loadFastaFile(tmpSeqFa)
    seqTbl.index()


def querySrcPsls(transMapSrcDbConn):
    srcAlignTbl = PslLite(transMapSrcDbConn, SourceDbTables.srcAlignsTbl)
    return srcAlignTbl.query("SELECT {} FROM {{table}};".format(PslLite.columnsNamesSql))
