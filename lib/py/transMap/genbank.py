from pycbio.hgdata import hgDb
import pipettor
from transMap.genomeDefs import AnnSetType
from transMap.transMapLite import TransMapSrcGene, getAccvSubselectClause


def valOrNone(val):
    return None if (val == "") or (val == "n/a") else val


def strOrNoneIfZero(val):
    return None if val == 0 else str(val)


class GenbankHgData(object):
    "obtain data UCSC browser data from various GenBank or RefSeq tables"

    pslSelectTmpl = """SELECT matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, """ \
                    """tBaseInsert, strand, concat(qName, ".", version), qSize, qStart, qEnd, tName, tSize, tStart, tEnd, """ \
                    """blockCount, blockSizes, qStarts, tStarts from {}, gbCdnaInfo """ \
                    """WHERE (qName = acc)"""

    def __init__(self, srcHgDb, annSetType):
        self.srcHgDb = srcHgDb
        self.annSetType = annSetType

    def __getPslTbl(self):
        if self.annSetType == AnnSetType.mrna:
            return "all_mrna"
        elif self.annSetType == AnnSetType.splicedEst:
            return "intronEst"
        elif self.annSetType == AnnSetType.refSeq:
            return "refSeqAli"

    def alignReader(self, limit=None):
        """get generator to return alignments with updated qNames,
        these are raw rows."""
        sql = self.pslSelectTmpl.format(self.__getPslTbl())
        if limit is not None:
            sql += " LIMIT {}".format(limit)
        hgsqlCmd = ("hgsql", "-Ne", sql, self.srcHgDb)
        sortCmd = ("sort", "-k", "14,14", "-k", "16,16n")
        uniqCmd = ("pslQueryUniq", "-p", "{}:".format(self.srcHgDb))

        with pipettor.Popen([hgsqlCmd, sortCmd, uniqCmd]) as fh:
            for line in fh:
                yield line.rstrip().split("\t")

    @staticmethod
    def __getTestSubsetClause(testAccvSubset):
        if testAccvSubset is None:
            return ""
        else:
            return "WHERE {}".format(getAccvSubselectClause("qName", testAccvSubset))

    def __metadataRowGen(self, sql, rowFactory, rowFilter=None):
        conn = hgDb.connect(self.srcHgDb, dictCursor=True)
        try:
            cur = conn.cursor()
            cur.execute(sql)
            for row in cur:
                if (rowFilter is None) or rowFilter(row):
                    yield rowFactory(row)
        finally:
            conn.close()

    def __refSeqToMetadata(self, row):
        geneType = "protein_coding" if row["accv"].startswith("NM_") else "non_coding"
        return TransMapSrcGene(srcId="{}:{}".format(self.srcHgDb, row["accv"]),
                               accv=row["accv"],
                               cds=valOrNone(row["cds_name"]),
                               geneId=strOrNoneIfZero(row["rl_locusLinkId"]),
                               geneName=valOrNone(row["rl_name"]),
                               geneType=geneType,
                               transcriptType=geneType)

    def __refSeqMetadataReader(self, testAccvSubset):
        sql = """SELECT concat(gb.acc, ".", gb.version) as accv, cds.name as cds_name, """ \
              """        rl.name as rl_name, rl.locusLinkId as rl_locusLinkId """ \
              """FROM hgFixed.gbCdnaInfo gb, hgFixed.cds, hgFixed.refLink rl  """ \
              """WHERE (gb.acc = rl.mrnaAcc) and (cds.id = gb.cds) and """ \
              """(gb.acc in (SELECT qName FROM refSeqAli {testAccvSubsetClause}));""".format(testAccvSubsetClause=self.__getTestSubsetClause(testAccvSubset))
        return self.__metadataRowGen(sql, self.__refSeqToMetadata)

    @staticmethod
    def __mrnaRowFilter(row):
        "only output rows with some data"
        return (row["cds_name"] != "n/a") or (row["gn_name"] != "n/a")

    def __mrnaToSrcMetadata(self, row):
        geneType = "protein_coding" if row["cds_name"] != "n/a" else "unknown"
        return TransMapSrcGene(srcId="{}:{}".format(self.srcHgDb, row["accv"]),
                               accv=row["accv"],
                               cds=valOrNone(row["cds_name"]),
                               geneId=None,
                               geneName=valOrNone(row["gn_name"]),
                               geneType=geneType,
                               transcriptType=geneType)

    def __refMRnaMetadataReader(self, testAccvSubset):
        sql = """select concat(gb.acc, ".", gb.version) as accv, cds.name as cds_name, gn.name as gn_name """ \
              """FROM hgFixed.gbCdnaInfo gb, hgFixed.cds, hgFixed.geneName gn """ \
              """WHERE (cds.id = gb.cds) and (gn.id = gb.geneName) and """ \
              """(gb.acc in (SELECT qName from all_mrna {testAccvSubsetClause}));""".format(testAccvSubsetClause=self.__getTestSubsetClause(testAccvSubset))
        return self.__metadataRowGen(sql, self.__mrnaToSrcMetadata, self.__mrnaRowFilter)

    def metadataReader(self, testAccvSubset=None):
        "reader for metadata for type; not valid for ESTs"
        if self.annSetType == AnnSetType.mrna:
            return self.__refMRnaMetadataReader(testAccvSubset)
        elif self.annSetType == AnnSetType.splicedEst:
            raise Exception("not for ESTs")
        elif self.annSetType == AnnSetType.refSeq:
            return self.__refSeqMetadataReader(testAccvSubset)
