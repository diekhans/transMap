import re
from pycbio.sys import fileOps
from .genomeDefs import CDnaType

class GenbankConf(object):
    """configuration of genbank for a genome or the defaults.
    handles on/off/defaulted"""
    __slots__ = ("db", "defaultConf", "on", "off")
    def __init__(self, db, defaultConf):
        self.db = db
        self.defaultConf = defaultConf
        self.on = set()
        self.off = set()

    def __str__(self):
        return self.db + "\ton:" + setOps.setJoin(self.on,",") + "\toff:" + setOps.setJoin(self.off,",")

    def getEnabled(self):
        return frozenset(self.on | self.defaultConf.on.difference(self.off))

class GenbankConfTbl(dict):
    "genbank config, indexed by database"

    def __init__(self):
        self.defaultConf = GenbankConf("default", None)
        self.__load()

    # parse clusterGenome line, used so we find entries that all all defaults
    parseGenomeRe = re.compile("^([A-Za-z0-9]+)\\.clusterGenome\s*=")

    # parse a cDNA line        1                2                  3                     4                  5
    parseCDnaRe = re.compile("^([A-Za-z0-9]+)\\.(genbank|refseq)\\.(mrna|est)\\.native\\.(load|align)\s*=\s*(yes|no)")

    def __load(self):
        "build dict of Conf objects"
        for line in fileOps.iterLines("/hive/data/outside/genbank/etc/genbank.conf"):
            self.__parseRow(line.strip())

    def __parseRow(self, line):
        geneomeMatch = self.parseGenomeRe.match(line)
        if geneomeMatch != None:
            self.__obtainGenbankConf(geneomeMatch.group(1))
        cdnaMatch = self.parseCDnaRe.match(line)
        if cdnaMatch != None:
            self.__processGenbankConf(cdnaMatch, line)

    def __obtainGenbankConf(self, db):
        "get a Conf entry, creating if it doesn't exist"
        if db == "default":
            return self.defaultConf
        else:
            if db not in self:
                self[db] = GenbankConf(db, self.defaultConf)
            return self[db]
            
    def __processGenbankConf(self, cdnaMatch, line):
        "add or update a GenbankConf object given a parsed line"
        db = cdnaMatch.group(1)
        cdna = None
        if cdnaMatch.group(2) == "refseq":
            cdna = CDnaType.refSeq
        elif cdnaMatch.group(2) == "genbank":
            if cdnaMatch.group(3) == "mrna":
                cdna = CDnaType.mrna
            elif cdnaMatch.group(3) == "est":
                cdna = CDnaType.splicedEst
        state = None
        if cdnaMatch.group(5) == "yes":
            state = True
        elif cdnaMatch.group(5) == "no":
            state = False
        if (cdna == None) or (state == None):
            raise Exception("can't parse genbank.conf line: "+line)
        conf = self.__obtainGenbankConf(db)
        if state:
            conf.on.add(cdna)
        else:
            conf.off.add(cdna)


