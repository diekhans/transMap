"""ExRun rule classes to obtain cDNA data for mappings.
"""
from pycbio.exrun import *
from transMap.GenomeDefs import CDnaTypes
from transMap import TransMap

def mkGbGetSeqsCmd(tm, srcDb, cdnaType, what, accFile=None):
    "what is argument to -get="
    cmd=["/cluster/data/genbank/bin/" + tm.arch + "/gbGetSeqs",
         "-db=" + srcDb.db,
         "-get=" + what,
         "-native", "-inclVersion",
         "-gbRoot=/cluster/data/genbank"]
    if accFile != None:
        cmd.append("-accFile="+accFile)
    cmd.extend((("refseq" if (cdnaType == CDnaTypes.refSeq) else "genbank"),
                ("est" if (cdnaType == CDnaTypes.splicedEst) else "mrna"),
                "stdout"))
    return cmd

class CDnaAligns(CmdRule):
    "obtain cDNA PSL alignments, making the alignment ids unique "
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        psl = tm.getSrcPsl(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType,
                              ("intronPsl" if (cdnaType == CDnaTypes.splicedEst) else "psl"))
        cmd2 = ("pslQueryUniq",)
        CmdRule.__init__(self, Cmd((cmd1, cmd2), stdout=psl))

class CDnaAlignStats(CmdRule):
    "obtain cDNA PSL alignments statistics"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        psl = tm.getSrcPsl(srcDb, cdnaType)
        stats = tm.getSrcAlnStats(srcDb, cdnaType)
        CmdRule.__init__(self, Cmd(("pslStats", IFileRef(psl), OFileRef(stats))))

class CDnaSeqIds(CmdRule):
    "obtain list of cDNA sequence ids (without unique suffix)"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        psl = tm.getSrcPsl(srcDb, cdnaType)
        ids = tm.getSrcSeqId(srcDb, cdnaType)
        CmdRule.__init__(self, Cmd((("bzcat", IFileRef(psl)),
                                    ("cut", "-f", "10"),
                                    ("sed", "-re", "s/-.+$//"),
                                    ("sort", "-u")),
                                   stdout=ids))

class CDnaSeqs(CmdRule):
    "obtain cDNA sequences"
    def __init__(self, tm, srcDb, cdnaType):
        self.tm = tm
        ids = tm.getSrcSeqId(srcDb, cdnaType)
        fa = tm.getSrcFa(srcDb, cdnaType)
        cmd = mkGbGetSeqsCmd(tm, srcDb, cdnaType, "seq", accFile=IFileRef(ids))
        CmdRule.__init__(self, Cmd(cmd, stdout=fa))

class CDnaCds(CmdRule):
    "obtain CDS specs"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        ids = tm.getSrcSeqId(srcDb, cdnaType)
        cds = tm.getSrcCds(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType, "ra")
        
        awkProg = """
        BEGIN {OFS = "\t"}
        /^acc / {acc = $2}
        /^cds / {cds = $2}
        /^$/ && (cds!="") {print acc,cds}
        /^$/ {acc=cds=""}
        """
        cmd2 = ("awk", awkProg)
        CmdRule.__init__(self, Cmd((cmd1, cmd2), stdout=cds))

class CDnaPolyA(CmdRule):
    "obtain poly information on cDNA sequences"
    # not currently used
    def __init__(self, tm, srcDb, cdnaType):
        self.tm = tm
        fa = tm.getSrcFa(srcDb, cdnaType)
        polya = tm.getSrcPoly(srcDb, cdnaType)
        cmd = ["faPolyASizes", IFileRef(fa), OFileRef(polya)]
        CmdRule.__init__(self, Cmd(cmd))


class UcscGenesAligns(CmdRule):
    "obtain PSLs and CDS files for UCSC genes"
    def __init__(self, tm, srcDb):
        self.tm = tm
        psl = tm.getSrcPsl(srcDb, CDnaTypes.ucscGenes)
        cds = tm.getSrcCds(srcDb, CDnaTypes.ucscGenes)
        
        cmd1 = ("hgsql", "-Ne",
                "select name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds from knownGene",
                srcDb.db)
        cmd2 = ("genePredToFakePsl", srcDb.db, "stdin", OFileRef(psl), OFileRef(cds))
        CmdRule.__init__(self, Cmd((cmd1, cmd2)))
    
class UcscGenesSeqs(CmdRule):
    def __init__(self, tm, srcDb):
        self.tm = tm
        fa = tm.getSrcFa(srcDb, CDnaTypes.ucscGenes)
        cmd1 = ("hgsql", "-Ne",
                """select concat(">",name,"\n",seq) from knownGeneMrna""",
                srcDb.db)
        CmdRule.__init__(self, Cmd((cmd1,), stdout=fa))

def createRules(tm, srcDb):
    "create rules to get data files for srcDb"
    for cdnaType in srcDb.cdnaTypes:
        if cdnaType == CDnaTypes.ucscGenes:
            tm.exrun.addRule(UcscGenesAligns(tm, srcDb))
            tm.exrun.addRule(UcscGenesSeqs(tm, srcDb))
        else:
            tm.exrun.addRule(CDnaAligns(tm, srcDb, cdnaType))
            tm.exrun.addRule(CDnaSeqs(tm, srcDb, cdnaType))

        tm.exrun.addRule(CDnaAlignStats(tm, srcDb, cdnaType))
        tm.exrun.addRule(CDnaSeqIds(tm, srcDb, cdnaType))
        if cdnaType in (CDnaTypes.refSeq, CDnaTypes.mrna):
            tm.exrun.addRule(CDnaCds(tm, srcDb, cdnaType))
