import re
from pycbio.sys import fileOps
from .genomeDefs import AnnSetType

# FIXME: drop unused
def mkGbGetSeqsCmd(srcDbName, annSetType, what, accFile=None):
    """create command to run gbGetSeqs.
    what is the argument to -get="""
    cmd=["/hive/data/outside/genbank/bin/x86_64/gbGetSeqs",
         "-db=" + srcDbName,
         "-get=" + what,
         "-native", "-inclVersion",
         "-gbRoot=/hive/data/outside/genbank"]
    if accFile != None:
        cmd += ["-accFile="+accFile]
    cmd += ["refseq" if (annSetType == AnnSetType.refSeq) else "genbank",
            "est" if (annSetType == AnnSetType.splicedEst) else "mrna",
            "/dev/stdout"]
    return cmd

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

    # parse clusterGenome line, used so we find entries that are all defaults
    parseGenomeRe = re.compile("^([A-Za-z0-9]+)\\.clusterGenome\s*=")

    # parse a cDNA line        1                2                  3                     4                  5
    parseAnnSetRe = re.compile("^([A-Za-z0-9]+)\\.(genbank|refseq)\\.(mrna|est)\\.native\\.(load|align)\s*=\s*(yes|no)")

    def __load(self):
        "build dict of Conf objects"
        for line in fileOps.iterLines("/hive/data/outside/genbank/etc/genbank.conf"):
            self.__parseRow(line.strip())

    def __parseRow(self, line):
        geneomeMatch = self.parseGenomeRe.match(line)
        if geneomeMatch != None:
            self.__obtainGenbankConf(geneomeMatch.group(1))
        annSetMatch = self.parseAnnSetRe.match(line)
        if annSetMatch != None:
            self.__processGenbankConf(annSetMatch, line)

    def __obtainGenbankConf(self, db):
        "get a Conf entry, creating if it doesn't exist"
        if db == "default":
            return self.defaultConf
        else:
            if db not in self:
                self[db] = GenbankConf(db, self.defaultConf)
            return self[db]
            
    def __processGenbankConf(self, annSetMatch, line):
        "add or update a GenbankConf object given a parsed line"
        db = annSetMatch.group(1)
        annSetType = None
        if annSetMatch.group(2) == "refseq":
            annSetType = AnnSetType.refSeq
        elif annSetMatch.group(2) == "genbank":
            if annSetMatch.group(3) == "mrna":
                annSetType = AnnSetType.mrna
            elif annSetMatch.group(3) == "est":
                annSetType = AnnSetType.splicedEst
        state = None
        if annSetMatch.group(5) == "yes":
            state = True
        elif annSetMatch.group(5) == "no":
            state = False
        if (annSetType == None) or (state == None):
            raise Exception("can't parse genbank.conf line: "+line)
        conf = self.__obtainGenbankConf(db)
        if state:
            conf.on.add(annSetType)
        else:
            conf.off.add(annSetType)


