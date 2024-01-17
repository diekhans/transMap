import os
import pipettor
from pycbio.hgdata import hgDb
from transMap import setSortLocale
from transMap.genomeData import AnnotationType
from transMap.srcData import SrcMetadata, getAccvSubselectClause
from pycbio.sys import fileOps
from pycbio.db import mysqlOps


def valOrNone(val):
    return None if (val == "") or (val == "n/a") else val


def strOrNoneIfZero(val):
    return None if val == 0 else str(val)


def haveRnaTrack(hgDbConn, hgDb):
    "check if all_mrna exists"
    return mysqlOps.haveTablesLike(hgDbConn, "all_mrna", hgDb)


def haveEstTrack(hgDbConn, hgDb):
    "check if intronEst exists"
    # might be split
    return mysqlOps.haveTablesLike(hgDbConn, "%intronEst", hgDb)


def haveRefseqTrack(hgDbConn, hgDb):
    "check if ncbiRefSeqPsl exists"
    return mysqlOps.haveTablesLike(hgDbConn, "ncbiRefSeqPsl", hgDb)


class GenbankHgData(object):
    "obtain data UCSC browser data from various GenBank or RefSeq tables"

    gbPslSelectTmpl = "SELECT matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, " \
        "tBaseInsert, strand, concat(qName, '.', version), qSize, qStart, qEnd, tName, tSize, tStart, tEnd, " \
        "blockCount, blockSizes, qStarts, tStarts FROM {}, hgFixed.gbCdnaInfo " \
        "WHERE (qName = acc)"

    # some refseq tables has non-accession entries (e.g. rna7), so just ignore those.
    rsPslSelectTmpl = "SELECT matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, " \
        "tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, " \
        "blockCount, blockSizes, qStarts, tStarts FROM {} WHERE qName REGEXP '^[A-Z0-9_]+\\\\.[0-9]+$'"

    def __init__(self, srcHgDb, annotationType):
        self.srcHgDb = srcHgDb
        self.annotationType = annotationType

    def _getIntronEstTables(self):
        hgDbConn = hgDb.connect(self.srcHgDb)
        try:
            return mysqlOps.getTablesLike(hgDbConn, "%intronEst", self.srcHgDb)
        finally:
            hgDbConn.close()

    def _getPslTbls(self):
        """get the list of tables for a source; intronEst might be split"""
        if self.annotationType == AnnotationType.rna:
            return ["all_mrna"]
        elif self.annotationType == AnnotationType.est:
            return self._getIntronEstTables()
        elif self.annotationType == AnnotationType.refseq:
            return ["ncbiRefSeqPsl"]

    def _buildPslSelect(self, table, testAccSubset=None):
        if self.annotationType == AnnotationType.refseq:
            sql = self.rsPslSelectTmpl.format(table)
        else:
            sql = self.gbPslSelectTmpl.format(table)
        if testAccSubset is not None:
            sql += " AND (qName in ({}))".format(",".join(['"{}"'.format(name) for name in testAccSubset]))
        return sql

    def _getPslCmdSingle(self, table, testAccSubset=None):
        hgsqlCmd = ("hgsql", "-Ne", self._buildPslSelect(table, testAccSubset), self.srcHgDb)
        return hgsqlCmd, None

    def _getPslCmdSplit(self, tables, testAccSubset=None):
        tmpPslFile = fileOps.tmpFileGet("transMap.splittable.", "psl")
        with open(tmpPslFile, "w") as pslFh:
            for table in tables:
                pipettor.run(("hgsql", "-Ne", self._buildPslSelect(table, testAccSubset), self.srcHgDb),
                             stdout=pslFh)
        return ["cat", tmpPslFile], tmpPslFile

    def alignReader(self, testAccSubset=None):
        """get generator to return alignments with updated qNames,
        these are raw rows."""
        setSortLocale()
        pslTables = self._getPslTbls()
        if len(pslTables) == 1:
            getAlignCmd, tmpPslFile = self._getPslCmdSingle(pslTables[0], testAccSubset)
        else:
            getAlignCmd, tmpPslFile = self._getPslCmdSplit(pslTables, testAccSubset)
        sortCmd = ("sort", "-k", "14,14", "-k", "16,16n")
        uniqCmd = ("pslQueryUniq", "--prefix", "{}:".format(self.srcHgDb))
        with pipettor.Popen([getAlignCmd, sortCmd, uniqCmd]) as fh:
            for line in fh:
                yield line.rstrip().split("\t")
        if tmpPslFile is not None:
            os.unlink(tmpPslFile)

    @staticmethod
    def _getTestSubsetClause(alignField, testAccvSubset):
        if testAccvSubset is None:
            return ""
        elif len(testAccvSubset) == 0:
            raise Exception("testAccvSubset is empty")
        else:
            return "AND {}".format(getAccvSubselectClause(alignField, testAccvSubset))

    def _metadataRowGen(self, sql, rowFactory):
        conn = hgDb.connect(self.srcHgDb, cursorclass=mysqlOps.DictCursor)
        try:
            cur = conn.cursor()
            cur.execute(sql)
            for row in cur:
                yield rowFactory(row)
        finally:
            conn.close()

    def _refseqToMetadata(self, row):
        geneType = "protein_coding" if row["accv"].startswith("NM_") else "non_coding"
        return SrcMetadata(srcId="{}:{}".format(self.srcHgDb, row["accv"]),
                           accv=row["accv"],
                           cds=valOrNone(row["cds_name"]),
                           geneId=strOrNoneIfZero(row["rl_locusLinkId"]),
                           geneName=valOrNone(row["rl_name"]),
                           geneType=geneType,
                           transcriptType=geneType)

    def _refseqMetadataReader(self, testAccvSubset):
        # using sub-select to get aligned sequence took forever, hence join and distinct
        sql = "SELECT DISTINCT psl.qName as accv, cds.cds as cds_name, " \
              "        rl.name as rl_name, rl.locusLinkId as rl_locusLinkId " \
              "FROM {srcHgDb}.ncbiRefSeqPsl psl " \
              "LEFT JOIN {srcHgDb}.ncbiRefSeqCds cds ON (cds.id = psl.qName) " \
              "LEFT JOIN {srcHgDb}.ncbiRefSeqLink rl ON (rl.mrnaAcc = psl.qName)  " \
              "WHERE TRUE {testAccvSubsetClause};" \
              .format(srcHgDb=self.srcHgDb,
                      testAccvSubsetClause=self._getTestSubsetClause("psl.qName", testAccvSubset))
        return self._metadataRowGen(sql, self._refseqToMetadata)

    def _rnaToSrcMetadata(self, row):
        geneType = "protein_coding" if row["cds_name"] != "n/a" else "unknown"
        return SrcMetadata(srcId="{}:{}".format(self.srcHgDb, row["accv"]),
                           accv=row["accv"],
                           cds=valOrNone(row["cds_name"]),
                           geneId=None,
                           geneName=valOrNone(row["gn_name"]),
                           geneType=geneType,
                           transcriptType=geneType)

    def _refRnaMetadataReader(self, testAccvSubset):
        # using sub-select to get aligned sequence took forever, hence join and distinct
        sql = """SELECT DISTINCT concat(gb.acc, ".", gb.version) as accv, cds.name as cds_name, gn.name as gn_name """ \
              """FROM hgFixed.gbCdnaInfo gb, hgFixed.cds, hgFixed.geneName gn, {srcHgDb}.all_mrna rna """ \
              """WHERE (cds.id = gb.cds) AND (gn.id = gb.geneName) AND (rna.qName =.gb.acc)""" \
              """{testAccvSubsetClause};""".format(srcHgDb=self.srcHgDb,
                                                   testAccvSubsetClause=self._getTestSubsetClause("rna.qName", testAccvSubset))
        return self._metadataRowGen(sql, self._rnaToSrcMetadata)

    def metadataReader(self, testAccvSubset=None):
        "reader for metadata for type; not valid for ESTs"
        if self.annotationType == AnnotationType.rna:
            return self._refRnaMetadataReader(testAccvSubset)
        elif self.annotationType == AnnotationType.est:
            raise Exception("not for ESTs")
        elif self.annotationType == AnnotationType.refseq:
            return self._refseqMetadataReader(testAccvSubset)
