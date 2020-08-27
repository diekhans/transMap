""""Database with genomes and chains."""
from collections import namedtuple, defaultdict
from pycbio.sys.symEnum import SymEnum
from pycbio.sys import fileOps
from pycbio.hgdata.hgSqlite import HgSqliteTable
from pycbio.db import sqliteOps
from transMap import hgDbNameParse


class GenomesDbTables(object):
    """tables names in a genomes sqlite3 database"""
    chainsTbl = "chains"
    genomeAsmsTbl = "genomeAsms"


class AnnotationType(SymEnum):
    "annotation type for a genome assembly"
    refseq = 0x01
    rna = 0x02
    est = 0x04
    ensembl = 0x08   # includes gencode


class AnnotationTypeSet(frozenset):
    "annotation types set for a genome assembly"

    def __str__(self):
        values = []
        for value in sorted(self):
            values.append(str(value))
        return ",".join(values)

    @staticmethod
    def fromStr(sval):
        if (sval is None) or (len(sval) == 0):
            return AnnotationTypeSet()
        else:
            return AnnotationTypeSet([AnnotationType(str(v)) for v in sval.split(',')])


class ChainType(SymEnum):
    "type of chains, which can either be existing or ones filtered on the fly."
    all = 1
    syn = 2
    rbest = 3


class Chains(namedtuple("Chains",
                        ("srcHgDb", "destHgDb", "chainType", "chainFile", "netFile"))):
    "Describes the chains used in a mapping to a destHgDb"
    __slots__ = ()

    @staticmethod
    def rowFactory(cur, row):
        return Chains(row[0], row[1], ChainType(str(row[2])), row[3], row[4])

    def toRow(self):
        # was sqlite3.register_adapter(ChainType, lambda chainType: str(chainType))
        return (self.srcHgDb, self.destHgDb, str(self.chainType), self.chainFile, self.netFile)


class ChainsSqliteTable(HgSqliteTable):
    """Interface todata on chains stored in the database"""
    _createSql = """CREATE TABLE {table} (
            srcHgDb text not null,
            destHgDb text not null,
            chainType text not null,
            chainFile text not null,
            netFile text not null);"""
    _insertSql = """INSERT INTO {table} (srcHgDb, destHgDb, chainType, chainFile, netFile) VALUES (?, ?, ?, ?, ?);"""
    _indexSql = ["""CREATE INDEX {table}_srcHgDb on {table} (srcHgDb);""",
                 """CREATE INDEX {table}_destHgDb on {table} (destHgDb);"""]

    columnNames = ("srcHgDb", "destHgDb", "chainType", "chainFile", "netFile")

    def __init__(self, conn, table, create=False):
        super().__init__(conn, table)
        if create:
            self._create(self._createSql)

    def index(self):
        """create index after loading"""
        self._index(self._indexSql)

    def loads(self, rows):
        """load rows into table"""
        self._inserts(self._insertSql, self.columnNames, [r.toRow() for r in rows])

    def queryByDbs(self, srcHgDb, destHgDb):
        sql = "SELECT {columns} FROM {table} WHERE (srcHgDb=?) AND (destHgDb=?)"
        return self.queryRows(sql, self.columnNames, Chains.rowFactory, srcHgDb, destHgDb)

    def queryByDestHgDb(self, destHgDb):
        sql = "SELECT {columns} FROM {table} WHERE (destHgDb=?)"
        return self.queryRows(sql, self.columnNames, Chains.rowFactory, destHgDb)

    def queryByDbsType(self, srcHgDb, destHgDb, chainType):
        "return Chain or None"
        sql = "SELECT {columns} FROM {table} WHERE (srcHgDb=?) AND (destHgDb=?) AND (chainType=?)"
        rows = list(self.queryRows(sql, self.columnNames, Chains.rowFactory, srcHgDb, destHgDb, str(chainType)))
        if len(rows) == 0:
            return None
        elif len(rows) == 1:
            return rows[0]
        else:
            raise Exception("too many rows returned")


class GenomeAsm(namedtuple("GenomeAsm",
                           ("hgDb", "clade", "commonName", "scientificName", "orgAbbrev", "annotationTypeSet"))):
    "Information about a genome assembly, stored in genome database "

    def __new__(cls, hgDb, clade, commonName, scientificName, orgAbbrev, annotationTypeSet):
        self = super().__new__(cls, hgDb, clade, commonName, scientificName, orgAbbrev, annotationTypeSet)
        self.dbNum = hgDbNameParse(hgDb)[1]
        return self

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(self.hgDb, self.clade, self.commonName, self.scientificName, str(self.annotationTypeSet))

    def isFinished(self):
        "is this genome considered finished"
        return self.commonName in ("Human", "Mouse")

    @staticmethod
    def rowFactory(cur, row):
        return GenomeAsm(row[0], row[1], row[2], row[3], row[4], AnnotationTypeSet.fromStr(row[5]))

    def toRow(self):
        # was sqlite3.register_adapter(AnnotationTypeSet, str)
        assert self.annotationTypeSet is not None, "None annotationTypeSet: {}".format(str(self))
        return (self.hgDb, self.clade, self.commonName, self.scientificName, self.orgAbbrev, str(self.annotationTypeSet))


class GenomeAsmsSqliteTable(HgSqliteTable):
    """Interface todata on chains stored in the database"""
    _createSql = """CREATE TABLE {table} (
            hgDb text not null,
            clade text not null,
            commonName text not null,
            scientificName text not null,
            orgAbbrev text not null,
            annotationTypeSet text);"""
    _insertSql = """INSERT INTO {table} (hgDb, clade, commonName, scientificName, orgAbbrev, annotationTypeSet) VALUES (?, ?, ?, ?, ?, ?);"""
    _indexSql = ["""CREATE INDEX {table}_hgDb on {table} (hgDb);""",
                 """CREATE INDEX {table}_commonName on {table} (commonName);""",
                 """CREATE INDEX {table}_annotationTypeSet on {table} (annotationTypeSet);"""]
    columnNames = ("hgDb", "clade", "commonName", "scientificName",
                   "orgAbbrev", "annotationTypeSet",)

    def __init__(self, conn, table, create=False):
        super().__init__(conn, table)
        if create:
            self._create(self._createSql)

    def index(self):
        """create index after loading"""
        self._index(self._indexSql)

    def loads(self, rows):
        """load rows into table"""
        self._inserts(self._insertSql, self.columnNames, [r.toRow() for r in rows])

    def queryByClades(self, clades):
        sql = "SELECT {{columns}} FROM {{table}} WHERE clade in ({})".format(",".join((len(clades) * ["?"])))
        return self.queryRows(sql, self.columnNames, GenomeAsm.rowFactory, *clades)

    def queryByHdDb(self, hgDb):
        sql = "SELECT {columns} FROM {table} WHERE hgDb = ?"
        return self.queryRows(sql, self.columnNames, GenomeAsm.rowFactory, hgDb)


class Organism(object):
    """Databases for a given organism. We use common name to tie assembly dbs together, since db prefix
    has changed between releases: Baboon papHam1 ->papAnu2
    """

    def __init__(self, commonName):
        self.commonName = commonName
        self.genomeAsms = []

    def add(self, genomeAsm):
        "add a genome assembly"
        self.genomeAsms.append(genomeAsm)

    def sort(self):
        "sort all entries, by reverse dbNum"
        self.genomeAsms.sort(key=lambda a: a.dbNum, reverse=True)


class Mapping(namedtuple("Mapping",
                         ("destHgDb", "srcHgDb", "annotationType", "chains"))):
    """Mapping to a destHgDb from a srcHgDb for an annotation type"""
    pass


def mappingsFilter(mappings, srcHgDbSubset, destHgDbSubset):
    "subset mappings based on command line srcHgDb/destHgDb restrictions"
    if srcHgDbSubset is not None:
        mappings = [m for m in mappings if m.srcHgDb in srcHgDbSubset]
    if destHgDbSubset is not None:
        mappings = [m for m in mappings if m.destHgDb in destHgDbSubset]
    return mappings


class DestDbMapping(namedtuple("DestDbMapping",
                               ("destHgDb", "annotationType"))):
    """Mapping to a destHgDb for an annotation type"""
    pass


class Genomes(object):
    """Object containing genome definitions, loaded from database for clades
    of interested. Also used to query relevant changes for mapping """
    def __init__(self, conf):
        self.conf = conf
        self.genomeAsms = dict()   # by hgDb
        self.orgs = dict()         # commmonName to Organism (-> GenomeDb)
        self.hgDbToOrg = dict()

        self.genomeDbConn = sqliteOps.connect(conf.genomeDb)
        self._loadGenomes()
        self._finish()
        self.chainsTbl = ChainsSqliteTable(self.genomeDbConn, GenomesDbTables.chainsTbl)

    def _loadGenomes(self):
        genomeAsmsTbl = GenomeAsmsSqliteTable(self.genomeDbConn, GenomesDbTables.genomeAsmsTbl)
        for genomeAsm in genomeAsmsTbl.queryByClades(self.conf.clades):
            self._addGenomeAsm(genomeAsm)

    def _addGenomeAsm(self, genomeAsm):
        "add a genome db"
        self.genomeAsms[genomeAsm.hgDb] = genomeAsm
        org = self._obtainOrganism(genomeAsm.commonName)
        org.add(genomeAsm)
        self.hgDbToOrg[genomeAsm.hgDb] = org

    def _obtainOrganism(self, commonName):
        org = self.orgs.get(commonName)
        if org is None:
            org = self.orgs[commonName] = Organism(commonName)
        return org

    def _finish(self):
        """Finish up, making all links"""
        for org in self.orgs.values():
            org.sort()

    def getAsmByName(self, hgDb):
        "lookup db name or None"
        return self.genomeAsms.get(hgDb)

    def findCandidateDestHgDbs(self):
        """Find list of candidate destHgDb.  This is the newest from the target clades,
        plus previous ones required in the config.  This may not actually be used
        if their are no chains.
        """
        candidateDestHgDbs = set()
        for org in self.orgs.values():
            candidateDestHgDbs.add(org.genomeAsms[0].hgDb)
        for hgDb in self.conf.requiredPreviousDestHgDbs:
            if self.getAsmByName(hgDb) is None:
                raise Exception("config requiredPreviousDestHgDbs {} is not a valid hgDb".format(hgDb))
            candidateDestHgDbs.add(hgDb)
        return candidateDestHgDbs

    def _getDestDbChains(self, destHgDb):
        return self.chainsTbl.queryByDestHgDb(destHgDb)

    def _getChainsBySrcOrg(self, destHgDb):
        "skip if srcHgDb not in table, these are outside of the specified clades"
        chainsBySrcOrg = defaultdict(list)
        for chain in self._getDestDbChains(destHgDb):
            if chain.srcHgDb in self.hgDbToOrg:
                chainsBySrcOrg[self.hgDbToOrg[chain.srcHgDb]].append(chain)
        return chainsBySrcOrg

    def _getPreferedChainTypes(self, destHgDb, srcHgDb):
        destClade = self.genomeAsms[destHgDb].clade
        srcClade = self.genomeAsms[srcHgDb].clade
        return self.conf.cladePreferedChains[frozenset([destClade, srcClade])]

    def _findChainMatchingPrefered(self, preferedChainType, srcChainsList):
        for srcChain in srcChainsList:
            if srcChain.chainType == preferedChainType:
                return srcChain
        return None

    def _getDestHgDbSrcHgPreferedChain(self, destHgDb, srcHgDb, srcChainsList):
        "get most appropriate chain for the clade differences"
        for preferedChainType in self._getPreferedChainTypes(destHgDb, srcHgDb):
            matchChain = self._findChainMatchingPrefered(preferedChainType, srcChainsList)
            if matchChain is not None:
                return matchChain
        raise Exception("can't find preferedChainType for {} ({}) to {} ({}) in {}"
                        .format(destHgDb, self.genomeAsms[destHgDb].clade,
                                srcHgDb, self.genomeAsms[srcHgDb].clade,
                                srcChainsList))

    def _getDestHgDbSrcHgDbMappings(self, destHgDb, srcHgDb, srcChainsList):
        "get mappings from dest to src with required chains"
        preferedChain = self._getDestHgDbSrcHgPreferedChain(destHgDb, srcHgDb, srcChainsList)
        # build for each annotation set in src
        srcAsm = self.getAsmByName(srcHgDb)
        for annotationType in srcAsm.annotationTypeSet:
            yield Mapping(destHgDb, srcHgDb, annotationType, preferedChain)

    def _getDestHgDbSrcOrgMappings(self, destHgDb, srcOrg, srcOrgChainsList):
        srcHgDb = srcOrg.genomeAsms[0].hgDb  # newest
        srcHgDbChains = [c for c in srcOrgChainsList if c.srcHgDb == srcHgDb]
        if len(srcHgDbChains) > 0:
            return list(self._getDestHgDbSrcHgDbMappings(destHgDb, srcHgDb, srcHgDbChains))
        else:
            return []

    def _getDestHgDbMappings(self, destHgDb):
        "get mappings to this destHgDb"
        destHgDbMappings = []
        chainsBySrcOrg = self._getChainsBySrcOrg(destHgDb)
        for srcOrg in chainsBySrcOrg.keys():
            destHgDbMappings.extend(self._getDestHgDbSrcOrgMappings(destHgDb, srcOrg, chainsBySrcOrg[srcOrg]))
        return destHgDbMappings

    def getCandidateMappings(self):
        "get sorted list of Mapping object for candidate database"
        candidateMappings = []
        for destHgDb in self.findCandidateDestHgDbs():
            candidateMappings.extend(self._getDestHgDbMappings(destHgDb))
        return sorted(candidateMappings, key=lambda m: m)

    @staticmethod
    def _mappingToDestDbMapping(mapping):
        return DestDbMapping(mapping.destHgDb, mapping.annotationType)

    def getCandidateDestDbMappings(self):
        "build sorted list of DestDbMapping for all candidate mappings"
        destDbMappings = set()
        for mapping in self.getCandidateMappings():
            destDbMappings.add(self._mappingToDestDbMapping(mapping))
        return sorted(destDbMappings, key=lambda m: m)

    def getCandidateSrcHgDbs(self):
        "build sorted list of srcHg for all candidate mappings"
        srcHgDbs = set()
        for mapping in self.getCandidateMappings():
            srcHgDbs.add(mapping.srcHgDb)
        return sorted(list(srcHgDbs))

    def dump(self, fh):
        "print for debugging purposes"
        fileOps.prRowv(fh, "GenomeDefs: genomes:")
        for org in self.orgs.values():
            fileOps.prRowv(fh, "", org.commonName)
            for db in org.dbs:
                fileOps.prRowv(fh, "", "", db)
