from pycbio.hgdata import hgDb
from pycbio.sys import fileOps
import pipettor


class EnsemblHgData(object):
    "obtain data UCSC browser data for GENCODE or Ensembl"

    gencodeCompTblBase = "wgEncodeGencodeComp"
    ensGeneTbl = "ensGene"

    def __init__(self, srcHgDb):
        self.srcHgDb = srcHgDb
        self.srcHgDbConn = hgDb.connect(srcHgDb, dictCursor=True)
        self.gencodeVersion = self.__getGencodeVersion()  # None if ensembl
        self.haveEnsembl = (self.__getTablesLike(self.ensGeneTbl) is not None)

    def haveData(self):
        return (self.gencodeVersion is not None) or self.haveEnsembl

    def __getTablesLike(self, pat):
        sql = "show tables like \"{}\"".format(pat)
        cur = self.srcHgDbConn.cursor()
        try:
            cur.execute(sql)
            return [row[0] for row in cur]
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
        return tbls[len(self.gencodeCompTblBase):]

    def __getGenePredTbl(self):
        if self.gencodeVersion is not None:
            return self.gencodeCompTblBase + self.gencodeVersion
        else:
            return self.ensGeneTbl

    def getTmpPslFile(self):
        genePredTbl = self.__getGenePredTbl()
        pslTmpFile = fileOps.tmpFileGet("transMap.{}.".format(genePredTbl), ".psl")
        toPslCmd = ("genePredToFakePsl", self.hgDb, genePredTbl, pslTmpFile, "/dev/null", "/dev/stdout")
        sortCmd = ("sort", "-k", "14,14", "-k", "16,16n")
        uniqCmd = ("pslQueryUniq", "-p", "{}:".format(self.srcHgDb))
        pipettor.run([toPslCmd, sortCmd, uniqCmd])
        return pslTmpFile

    def getCdsSpecs(self):
        genePredTbl = self.__getGenePredTbl()
        with pipettor.Popen(["genePredToFakePsl", self.hgDb, genePredTbl, "/dev/null", "/dev/stdout"]) as cdsFh:
            return [line.split('\t') for line in fileOps.iterLines(cdsFh)]
