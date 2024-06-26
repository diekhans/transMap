import re
from pycbio.sys import fileOps, setOps
from .genomeData import AnnotationType


class GenbankDbConf(object):
    """configuration of genbank for a genome or the defaults.
    handles on/off/defaulted"""
    __slots__ = ("hgDb", "defaultConf", "_on", "_off")

    def __init__(self, hgDb, defaultConf):
        self.hgDb = hgDb
        self.defaultConf = defaultConf
        # sets of annotation types turned on/off, doesn't include defaults!!
        self._on = set()
        self._off = set()

    def addOn(self, annotationType):
        self._on.add(annotationType)

    def addOff(self, annotationType):
        self._off.add(annotationType)

    def __str__(self):
        return "{}\ton:{}\toff:{}".format(self.hgDb, setOps.setJoin(self._on, ","),
                                          setOps.setJoin(self._off, ","))

    def getEnabled(self):
        return frozenset(self._on | self.defaultConf._on.difference(self._off))


class GenbankConf(dict):
    "genbank config, indexed by database"

    stdConfRaFile = "/hive/data/outside/genbank/etc/genbank.conf"

    def __init__(self, confRaFile=stdConfRaFile):
        self.defaultConf = GenbankDbConf("default", None)
        self._load(confRaFile)

    # parse serverGenome line, used so we find entries that are all defaults
    parseGenomeRe = re.compile("^([A-Za-z0-9]+)\\.serverGenome\\s*=")

    # parse a cDNA line        1                2                  3                     4                  5
    parseAnnSetRe = re.compile("^([A-Za-z0-9]+)\\.(genbank|refseq)\\.(mrna|est)\\.native\\.(load|align)\\s*=\\s*(yes|no)")

    def _load(self, confRaFile):
        "build dict of Conf objects"
        for line in fileOps.iterLines(confRaFile):
            self._parseRow(line.strip())

    def _parseRow(self, line):
        genomeMatch = self.parseGenomeRe.match(line)
        if genomeMatch is not None:
            self._obtainGenbankDbConf(genomeMatch.group(1))
        annSetMatch = self.parseAnnSetRe.match(line)
        if annSetMatch is not None:
            self._processGenbankDbConf(annSetMatch, line)

    def _obtainGenbankDbConf(self, hgDb):
        "get a Conf entry, creating if it doesn't exist"
        if hgDb == "default":
            return self.defaultConf
        else:
            if hgDb not in self:
                self[hgDb] = GenbankDbConf(hgDb, self.defaultConf)
            return self[hgDb]

    def _processGenbankDbConf(self, annSetMatch, line):
        "add or update a GenbankDbConf object given a parsed line"
        hgDb = annSetMatch.group(1)
        annotationType = None
        if annSetMatch.group(2) == "refseq":
            annotationType = AnnotationType.refseq
        elif annSetMatch.group(2) == "genbank":
            if annSetMatch.group(3) == "mrna":
                annotationType = AnnotationType.rna
            elif annSetMatch.group(3) == "est":
                annotationType = AnnotationType.est
        state = None
        if annSetMatch.group(5) == "yes":
            state = True
        elif annSetMatch.group(5) == "no":
            state = False
        if (annotationType is None) or (state is None):
            raise Exception("can't parse genbank.conf line: {}".format(line))
        conf = self._obtainGenbankDbConf(hgDb)
        if state:
            conf.addOn(annotationType)
        else:
            conf.addOff(annotationType)
