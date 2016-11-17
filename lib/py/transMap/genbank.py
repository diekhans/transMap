from pycbio.hgdata import hgDb
import pipettor
from transMap.genomeData import AnnotationType
from transMap.srcData import SrcMetadata, getAccvSubselectClause


def valOrNone(val):
    return None if (val == "") or (val == "n/a") else val


def strOrNoneIfZero(val):
    return None if val == 0 else str(val)


class GenbankHgData(object):
    "obtain data UCSC browser data from various GenBank or RefSeq tables"

    pslSelectTmpl = """SELECT matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, """ \
                    """tBaseInsert, strand, concat(qName, ".", version), qSize, qStart, qEnd, tName, tSize, tStart, tEnd, """ \
                    """blockCount, blockSizes, qStarts, tStarts from {}, hgFixed.gbCdnaInfo """ \
                    """WHERE (qName = acc)"""

    def __init__(self, srcHgDb, annotationType):
        self.srcHgDb = srcHgDb
        self.annotationType = annotationType

    def __getPslTbl(self):
        if self.annotationType == AnnotationType.rna:
            return "all_mrna"
        elif self.annotationType == AnnotationType.est:
            return "intronEst"
        elif self.annotationType == AnnotationType.refseq:
            return "refSeqAli"

    def alignReader(self, testAccSubset=None):
        """get generator to return alignments with updated qNames,
        these are raw rows."""
        sql = self.pslSelectTmpl.format(self.__getPslTbl())
        if testAccSubset is not None:
            sql += " AND (qName in ({}))".format(",".join(['"{}"'.format(name) for name in testAccSubset]))
        hgsqlCmd = ("hgsql", "-Ne", sql, self.srcHgDb)
        sortCmd = ("sort", "-k", "14,14", "-k", "16,16n")
        uniqCmd = ("pslQueryUniq", "-p", "{}:".format(self.srcHgDb))
        with pipettor.Popen([hgsqlCmd, sortCmd, uniqCmd]) as fh:
            for line in fh:
                yield line.rstrip().split("\t")

    @staticmethod
    def __getTestSubsetClause(alignField, testAccvSubset):
        if testAccvSubset is None:
            return ""
        else:
            return "AND {}".format(getAccvSubselectClause(alignField, testAccvSubset))

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

    def __refseqToMetadata(self, row):
        geneType = "protein_coding" if row["accv"].startswith("NM_") else "non_coding"
        return SrcMetadata(srcId="{}:{}".format(self.srcHgDb, row["accv"]),
                           accv=row["accv"],
                           cds=valOrNone(row["cds_name"]),
                           geneId=strOrNoneIfZero(row["rl_locusLinkId"]),
                           geneName=valOrNone(row["rl_name"]),
                           geneType=geneType,
                           transcriptType=geneType)

    def __refseqMetadataReader(self, testAccvSubset):
        # using sub-select to get aligned sequence took forever, hence join and distinct
        sql = """SELECT DISTINCT concat(gb.acc, ".", gb.version) as accv, cds.name as cds_name, """ \
              """        rl.name as rl_name, rl.locusLinkId as rl_locusLinkId """ \
              """FROM hgFixed.gbCdnaInfo gb, hgFixed.cds, hgFixed.refLink rl, {srcHgDb}.refSeqAli ali """ \
              """WHERE (gb.acc = rl.mrnaAcc) AND (cds.id = gb.cds) AND (ali.qName = gb.acc) """ \
              """{testAccvSubsetClause};""".format(srcHgDb=self.srcHgDb,
                                                   testAccvSubsetClause=self.__getTestSubsetClause("ali.qName",
                                                                                                   testAccvSubset))
        return self.__metadataRowGen(sql, self.__refseqToMetadata)

    @staticmethod
    def __rnaRowFilter(row):
        "only output rows with some data"
        return (row["cds_name"] != "n/a") or (row["gn_name"] != "n/a")

    def __rnaToSrcMetadata(self, row):
        geneType = "protein_coding" if row["cds_name"] != "n/a" else "unknown"
        return SrcMetadata(srcId="{}:{}".format(self.srcHgDb, row["accv"]),
                           accv=row["accv"],
                           cds=valOrNone(row["cds_name"]),
                           geneId=None,
                           geneName=valOrNone(row["gn_name"]),
                           geneType=geneType,
                           transcriptType=geneType)

    def __refRnaMetadataReader(self, testAccvSubset):
        # using sub-select to get aligned sequence took forever, hence join and distinct
        sql = """SELECT DISTINCT concat(gb.acc, ".", gb.version) as accv, cds.name as cds_name, gn.name as gn_name """ \
              """FROM hgFixed.gbCdnaInfo gb, hgFixed.cds, hgFixed.geneName gn, {srcHgDb}.all_mrna rna """ \
              """WHERE (cds.id = gb.cds) AND (gn.id = gb.geneName) AND (rna.qName =.gb.acc)""" \
              """{testAccvSubsetClause};""".format(srcHgDb=self.srcHgDb,
                                                   testAccvSubsetClause=self.__getTestSubsetClause("rna.qName",
                                                                                                   testAccvSubset))
        return self.__metadataRowGen(sql, self.__rnaToSrcMetadata, self.__rnaRowFilter)

    def metadataReader(self, testAccvSubset=None):
        "reader for metadata for type; not valid for ESTs"
        if self.annotationType == AnnotationType.rna:
            return self.__refRnaMetadataReader(testAccvSubset)
        elif self.annotationType == AnnotationType.est:
            raise Exception("not for ESTs")
        elif self.annotationType == AnnotationType.refseq:
            return self.__refseqMetadataReader(testAccvSubset)
