from pycbio.hgdata import hgDb
from pycbio.sys import fileOps
from transMap.genomeData import AnnotationType
from transMap.srcData import SrcMetadata, getAccvSubselectClause
import pipettor

# FIXME:  filter by biotype to removed small RNAs


def valOrNone(val):
    return None if (val == "") else val


class EnsemblHgData(object):
    "obtain data UCSC browser data for GENCODE or Ensembl"

    gencodeCompTblBase = "wgEncodeGencodeComp"
    ensGeneTbl = "ensGene"

    def __init__(self, srcHgDb, annotationType):
        self.srcHgDb = srcHgDb
        self.annotationType = annotationType
        self.srcHgDbConn = hgDb.connect(srcHgDb, dictCursor=True)
        self.gencodeVersion = self.__getGencodeVersion() if annotationType == AnnotationType.gencode else None

    def __getTablesLike(self, pat):
        sql = "show tables like \"{}\"".format(pat)
        key = "Tables_in_{} ({})".format(self.srcHgDb, pat)
        cur = self.srcHgDbConn.cursor()
        try:
            cur.execute(sql)
            return [row[key] for row in cur]
        finally:
            cur.close()

    def __getGencodeVersion(self):
        "return latest gencode version, parse from table names, or None"
        tbls = self.__getTablesLike(self.gencodeCompTblBase + "%")
        if len(tbls) == 0:
            return None
        # get most current version number by sorting by reverse string length,
        # then reverse value will get the latest
        tbls.sort(key=lambda t: (len(t), t))
        return tbls[0][len(self.gencodeCompTblBase):]

    def __getGenePredTbl(self):
        if self.gencodeVersion is not None:
            return self.gencodeCompTblBase + self.gencodeVersion
        else:
            return self.ensGeneTbl

    def alignReader(self, limit=None):
        """get generator to return alignments with updated qNames,
        these are raw rows."""
        toPslCmd = ("genePredToFakePsl", self.srcHgDb, self.__getGenePredTbl(), "/dev/stdout", "/dev/null")
        sortCmd = ("sort", "-k", "14,14", "-k", "16,16n")
        uniqCmd = ("pslQueryUniq", "-p", "{}:".format(self.srcHgDb))
        cmd = [toPslCmd]
        if limit is not None:
            cmd.append(["tawk", "NR<={}".format(limit)])
        cmd.extend([sortCmd, uniqCmd])
        with pipettor.Popen(cmd, "r") as fh:
            for line in fh:
                yield line.rstrip().split("\t")

    def __processCdsSpecLine(self, line, testAccvSubset, cdsSpecs):
        accv, cds = line.split("\t")
        if (testAccvSubset is None) or (accv in testAccvSubset):
            if accv in cdsSpecs:
                if cds != cdsSpecs[accv]:
                    raise Exception("multiple occurrences of CDS for {} first was {}, second is {}".format(accv, cdsSpecs[accv], cds))
            else:
                cdsSpecs[accv] = cds

    def __getCdsSpecs(self, testAccvSubset):
        toCdsCmd = ["genePredToFakePsl", self.srcHgDb, self.__getGenePredTbl(), "/dev/null", "/dev/stdout"]
        with pipettor.Popen([toCdsCmd], "r") as cdsFh:
            cdsSpecs = {}
            for line in fileOps.iterLines(cdsFh):
                self.__processCdsSpecLine(line, testAccvSubset, cdsSpecs)
            return cdsSpecs

    def __metadataRowGen(self, sql, cdsSpecs, rowFactory):
        conn = hgDb.connect(self.srcHgDb, dictCursor=True)
        try:
            cur = conn.cursor()
            cur.execute(sql)
            for row in cur:
                yield rowFactory(cdsSpecs, row)
        finally:
            conn.close()

    def __toMetadata(self, cdsSpecs, row):
        cds = cdsSpecs.get(row["transcriptId"])
        return SrcMetadata(srcId="{}:{}".format(self.srcHgDb, row["transcriptId"]),
                           accv=row["transcriptId"],
                           cds=cds,
                           geneId=row["geneId"],
                           geneName=row["geneName"],
                           geneType=row["geneType"],
                           transcriptType=row["transcriptType"])

    def __gencodeMetadataReader(self, cdsSpecs, testAccvSubset):
        if testAccvSubset is not None:
            transcriptIdSelect = getAccvSubselectClause("transcriptId", testAccvSubset)
        else:
            transcriptIdSelect = """transcriptId in (SELECT name from wgEncodeGencodeComp{ver})""".format(ver=self.gencodeVersion)
        sql = """SELECT transcriptId, geneId, geneName, geneType, transcriptType """ \
              """FROM wgEncodeGencodeAttrs{ver} """ \
              """WHERE {transcriptIdSelect}""".format(ver=self.gencodeVersion, transcriptIdSelect=transcriptIdSelect)
        return self.__metadataRowGen(sql, cdsSpecs, self.__toMetadata)

    def __ensemblMetadataReader(self, cdsSpecs, testAccvSubset):
        if testAccvSubset is not None:
            transcriptIdSelect = " AND " + getAccvSubselectClause("eg.name", testAccvSubset)
        else:
            transcriptIdSelect = ""
        sql = """SELECT eg.name as transcriptId, eg.name2 as geneId, gn.value as geneName, es.source as geneType, es.source as transcriptType """ \
              """FROM ensGene eg LEFT JOIN ensemblToGeneName gn ON gn.name=eg.name, ensemblSource es """ \
              """WHERE (es.name = eg.name) {}""".format(transcriptIdSelect)
        print len(testAccvSubset), testAccvSubset
        print sql
        return self.__metadataRowGen(sql, cdsSpecs, self.__toMetadata)

    def metadataReader(self, testAccvSubset=None):
        "reader for metadata"
        cdsSpecs = self.__getCdsSpecs(testAccvSubset)
        if self.gencodeVersion is not None:
            return self.__gencodeMetadataReader(cdsSpecs, testAccvSubset)
        else:
            return self.__ensemblMetadataReader(cdsSpecs, testAccvSubset)
