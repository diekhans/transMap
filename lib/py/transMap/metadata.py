from pycbio.tsv import TabFileReader
from pycbio.sys import fileOps, dbOps
from pycbio.hgdata import hgDb

class TranscriptMetadata(object):
    def __init__(self, id):
        self.id = id
        self.cds = self.gene = self.cat = ""

    def haveData(self):
        return (self.cds != "") or (self.gene != "") or (self.cat != "")

class TranscriptMetadataTbl(dict):
    def loadCdsFile(self, cdsFile):
        for row in TabFileReader(cdsFile, hashAreComments=True):
            self.obtain(row[0]).cds = row[1]

    def __iterRows(self, conn, sql):
        for row in dbOps.queryConnection(conn, sql):
            yield row

    def loadGenesDb(self, conn, sql):
        for row in self.__iterRows(conn, sql):
            self.obtain(row[0]).gene = row[1]

    def loadCatDb(self, conn, sql):
        for row in self.__iterRows(conn, sql):
            self.obtain(row[0]).cat = row[1]

    def obtain(self, id):
        meta = self.get(id)
        if meta == None:
            meta = self[id] = TranscriptMetadata(id)
        return meta

    def write(self, outFh):
        ids = list(self.iterkeys())
        ids.sort()
        for id in ids:
            meta = self[id]
            if meta.haveData():
                fileOps.prRowv(outFh, id, meta.cds, meta.gene, meta.cat)
