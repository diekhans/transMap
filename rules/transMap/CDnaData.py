"""ExRun rule classes to obtain cDNA data for mappings.
"""
from transMap.GenomeDefs import CDnaTypes
from transMap import TransMap

def mkGbGetSeqsCmd(tm, srcDb, cdnaType, what):
    "what is argument to -get="
    cmd=("/cluster/data/genbank/bin/" + tm.arch + "/gbGetSeqs",
         "-db=" + srcDb.db,
         "-get=" + what,
         "-native", "-inclVersion",
         "-gbRoot=/cluster/data/genbank",
         ("refseq" if (cdnaType == CDnaTypes.refSeq) else "genbank"),
         ("est" if (cdnaType == CDnaTypes.splicedEst) else "mrna"),
         "stdout")
    return cmd

class CDnaAligns(CmdRule):
    "obtain cDNA PSL alignments, making the alignment ids unique "
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        self.tm = tm
        psl = tm.getSrcPsl(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType,
                              ("introlPsl" if (cdnaType == CDnaTypes.splicedEst) else "psl"))
        cmd2 = ("bin/pslQueryUniq",)
        cmd3 = ("bzip2", "-c")
        CmdRule.__init__(self, psl.realPath, (cmd1, cmd2, cmd3), stdout=psl)

class CDnaAlignStats(CmdRule):
    "obtain cDNA PSL alignments statistics"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        self.tm = tm
        psl = tm.getSrcPsl(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType,
                              ("introlPsl" if (cdnaType == CDnaTypes.splicedEst) else "psl"))
        cmd2 = ("bin/pslQueryUniq",)
        cmd3 = ("bzip2", "-c")
        CmdRule.__init__(self, psl.realPath, (cmd1, cmd2, cmd3), stdout=psl)

class CDnaSeqIds(CmdRule):
    "obtain list of cDNA sequence ids (without unique suffix) alignments"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        self.tm = tm
        psl = tm.getSrcPsl(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType,
                              ("introlPsl" if (cdnaType == CDnaTypes.splicedEst) else "psl"))
        cmd2 = ("bin/pslQueryUniq",)
        cmd3 = ("bzip2", "-c")
        CmdRule.__init__(self, psl.realPath, (cmd1, cmd2, cmd3), stdout=psl)

class CDnaSeqs(CmdRule):
    "obtain cDNA sequences"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        self.tm = tm
        fa = tm.getSrcFa(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType,
                              ("introlPsl" if (cdnaType == CDnaTypes.splicedEst) else "psl"))
        cmd3 = ("bzip2", "-c")
        CmdRule.__init__(self, psl.realPath, (cmd1, cmd2, cmd3), stdout=psl)

