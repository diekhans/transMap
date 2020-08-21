import re
from pycbio.hgdata import hgDb
from pycbio.sys import fileOps
from pycbio.db import mysqlOps
from transMap import setSortLocale
from transMap.srcData import SrcMetadata, getAccvSubselectClause
import pipettor

# FIXME:  filter by biotype to removed small RNAs

gencodeCompTblBase = "wgEncodeGencodeComp"
ensGeneTbl = "ensGene"


def valOrNone(val):
    return None if (val == "") else val


def isEnsemblProjectionBuild(hgDbConn, srcHgDb):
    """does this ensembl appear to be a projection (full or mixed) gene build?
    this is determined by there being at least one transcript that is mapped multiple times"""
    sql = "SELECT count(*) cnt FROM (SELECT name, count(*) AS cnt FROM {srcHgDb}.{ensGeneTbl} GROUP BY name) AS counts WHERE counts.cnt > 2".format(srcHgDb=srcHgDb, ensGeneTbl=ensGeneTbl)
    rows = list(mysqlOps.query(hgDbConn, sql))
    return (rows[0]["cnt"] > 0)


def haveEnsemblFull(hgDbConn, srcHgDb):
    """does this database have an Ensembl full or mixed build? (not projection"""
    return (mysqlOps.haveTablesLike(hgDbConn, ensGeneTbl, db=srcHgDb)
            and not isEnsemblProjectionBuild(hgDbConn, srcHgDb))


def haveGencode(hgDbConn, srcHgDb):
    return mysqlOps.haveTablesLike(hgDbConn, "{}%".format(gencodeCompTblBase), db=srcHgDb)


class EnsemblHgData(object):
    """Obtain data UCSC browser data for GENCODE or Ensembl.
    Optionally force a specific GENCODE version to avoid pre-releases"""

    def __init__(self, srcHgDb, annotationType, gencodeVersion):
        self.srcHgDb = srcHgDb
        self.annotationType = annotationType
        self.srcHgDbConn = hgDb.connect(srcHgDb, dictCursor=True)
        if gencodeVersion is not None:
            self.gencodeVersion = gencodeVersion
        else:
            self.gencodeVersion = self._findGencodeVersion()  # None if not GENCODE

    def _parseGenbankVersion(self, table):
        """parse numeric version from table name, avoiding V24lift37 like names"""
        verRe = "^{}(VM?[0-9]+)$".format(gencodeCompTblBase)
        m = re.match(verRe, table)
        if m is None:
            return None
        else:
            return m.group(1)

    def _getGencodeCompTables(self):
        tbls = []
        for tbl in mysqlOps.getTablesLike(self.srcHgDbConn, "{}%".format(gencodeCompTblBase)):
            # skip lift tables
            if self._parseGenbankVersion(tbl) is not None:
                tbls.append(tbl)
        return tbls

    def _findGencodeVersion(self):
        "return latest gencode version, parse from table names, or None"
        tbls = self._getGencodeCompTables()
        if len(tbls) == 0:
            return None
        # get most current version number by sorting by reverse string length,
        # then reverse value will get the latest
        tbls.sort(key=lambda t: (len(t), t), reverse=True)
        return tbls[0][len(gencodeCompTblBase):]

    def _getGenePredTbl(self):
        if self.gencodeVersion is not None:
            return gencodeCompTblBase + self.gencodeVersion
        else:
            return ensGeneTbl

    def _getGenePredSql(self, testAccvSubset):
        sql = "SELECT name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, "\
              "exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames " \
              "FROM {}".format(self._getGenePredTbl())
        if testAccvSubset is not None:
            if len(testAccvSubset) == 0:
                raise Exception("empty testAccvSubset provided for {}".format(self.srcHgDb))
            sql += " WHERE {}".format(getAccvSubselectClause("name", testAccvSubset))
        return sql

    def _getGenePredToPslCmds(self, pslOut, cdsOut, testAccvSubset):
        getGenePredCmd = ("hgsql", "-Ne", self._getGenePredSql(testAccvSubset), self.srcHgDb)
        toPslCmd = ("genePredToFakePsl", self.srcHgDb, "/dev/stdin", pslOut, cdsOut)
        return [getGenePredCmd, toPslCmd]

    def alignReader(self, testAccvSubset):
        """get generator to return alignments with updated qNames,
        these are raw rows."""
        setSortLocale()
        getGenePredCmd, toPslCmd = self._getGenePredToPslCmds("/dev/stdout", "/dev/null", testAccvSubset)
        sortCmd = ("sort", "-k", "14,14", "-k", "16,16n")
        uniqCmd = ("pslQueryUniq", "--prefix", "{}:".format(self.srcHgDb))
        with pipettor.Popen([getGenePredCmd, toPslCmd, sortCmd, uniqCmd], "r") as fh:
            for line in fh:
                yield line.rstrip().split("\t")

    def _processCdsSpecLine(self, line, testAccvSubset, cdsSpecs):
        accv, cds = line.split("\t")
        if (testAccvSubset is None) or (accv in testAccvSubset):
            if accv in cdsSpecs:
                if cds != cdsSpecs[accv]:
                    raise Exception("multiple occurrences of CDS for {} first was {}, second is {}".format(accv, cdsSpecs[accv], cds))
            else:
                cdsSpecs[accv] = cds

    def _getCdsSpecs(self, testAccvSubset):
        getGenePredCmd, toCdsCmd = self._getGenePredToPslCmds("/dev/null", "/dev/stdout", testAccvSubset)
        with pipettor.Popen([getGenePredCmd, toCdsCmd], "r") as cdsFh:
            cdsSpecs = {}
            for line in fileOps.iterLines(cdsFh):
                self._processCdsSpecLine(line, testAccvSubset, cdsSpecs)
            return cdsSpecs

    def _metadataRowGen(self, sql, cdsSpecs, rowFactory):
        conn = hgDb.connect(self.srcHgDb, dictCursor=True)
        try:
            cur = conn.cursor()
            cur.execute(sql)
            for row in cur:
                yield rowFactory(cdsSpecs, row)
        finally:
            conn.close()

    def _toMetadata(self, cdsSpecs, row):
        cds = cdsSpecs.get(row["transcriptId"])
        return SrcMetadata(srcId="{}:{}".format(self.srcHgDb, row["transcriptId"]),
                           accv=row["transcriptId"],
                           cds=cds,
                           geneId=row["geneId"],
                           geneName=row["geneName"],
                           geneType=row["geneType"],
                           transcriptType=row["transcriptType"])

    def _gencodeMetadataReader(self, cdsSpecs, testAccvSubset):
        if testAccvSubset is not None:
            transcriptIdSelect = getAccvSubselectClause("transcriptId", testAccvSubset)
        else:
            transcriptIdSelect = """transcriptId in (SELECT name from wgEncodeGencodeComp{ver})""".format(ver=self.gencodeVersion)
        sql = """SELECT transcriptId, geneId, geneName, geneType, transcriptType """ \
              """FROM wgEncodeGencodeAttrs{ver} """ \
              """WHERE {transcriptIdSelect}""".format(ver=self.gencodeVersion, transcriptIdSelect=transcriptIdSelect)
        return self._metadataRowGen(sql, cdsSpecs, self._toMetadata)

    def _ensemblMetadataReader(self, cdsSpecs, testAccvSubset):
        if testAccvSubset is not None:
            transcriptIdSelect = " AND " + getAccvSubselectClause("eg.name", testAccvSubset)
        else:
            transcriptIdSelect = ""
        sql = """SELECT eg.name as transcriptId, eg.name2 as geneId, gn.value as geneName, es.source as geneType, es.source as transcriptType """ \
              """FROM ensGene eg LEFT JOIN ensemblToGeneName gn ON gn.name=eg.name, ensemblSource es """ \
              """WHERE (es.name = eg.name) {}""".format(transcriptIdSelect)
        return self._metadataRowGen(sql, cdsSpecs, self._toMetadata)

    def metadataReader(self, testAccvSubset=None):
        "reader for metadata"
        cdsSpecs = self._getCdsSpecs(testAccvSubset)
        if self.gencodeVersion is not None:
            return self._gencodeMetadataReader(cdsSpecs, testAccvSubset)
        else:
            return self._ensemblMetadataReader(cdsSpecs, testAccvSubset)
