"""ExRun rule classes to obtain cDNA data for mappings.
"""
import re
from pycbio.exrun import *
from transMap.genomeDefs import CDnaType
from transMap import transMap
from pycbio.sys import fileOps, procOps

# N.B. /dev/stdout is used instead of stdout, as it appears related to fasta corruption problem

def mkGbGetSeqsCmd(tm, srcDb, cdnaType, what, accFile=None):
    "what is argument to -get="
    cmd=["/hive/data/outside/genbank/bin/" + tm.arch + "/gbGetSeqs",
         "-db=" + srcDb.name,
         "-get=" + what,
         "-native", "-inclVersion",
         "-gbRoot=/hive/data/outside/genbank"]
    if accFile != None:
        cmd.append("-accFile="+accFile)
    cmd.extend((("refseq" if (cdnaType == CDnaType.refSeq) else "genbank"),
                ("est" if (cdnaType == CDnaType.splicedEst) else "mrna"),
                "/dev/stdout"))
    return cmd

# FIXME: column stuff is duplicated here and in commands
def makeColumnWithSrcDbConcat(db, column):
    return 'concat("' + db + ':", ' + column + ')'

class CDnaAligns(CmdRule):
    "obtain cDNA PSL alignments, making the alignment ids unique "
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        psl = tm.getSrcPsl(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType,
                              ("intronPsl" if (cdnaType == CDnaType.splicedEst) else "psl"))
        cmd2 = ("pslQueryUniq", "-p", srcDb.name+":")
        CmdRule.__init__(self, Cmd((cmd1, cmd2), stdout=psl))

class CDnaAlignStats(CmdRule):
    "obtain cDNA PSL alignments statistics"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        psl = tm.getSrcPsl(srcDb, cdnaType)
        stats = tm.getSrcAlnStats(srcDb, cdnaType)
        CmdRule.__init__(self, Cmd(("pslStats", FileIn(psl), FileOut(stats))))

class CDnaSeqIds(CmdRule):
    "obtain list of cDNA sequence ids (without srcDb prefix or unique suffix)"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        psl = tm.getSrcPsl(srcDb, cdnaType)
        ids = tm.getSrcSeqId(srcDb, cdnaType)
        CmdRule.__init__(self, Cmd((("cut", "-f", "10", FileIn(psl)),
                                    ("sed", "-r", "-e", "s/-.+$//",  "-e", "s/^" + srcDb.name+ "://"),
                                    ("sort", "-u")),
                                   stdout=ids))

class CDnaSeqs(CmdRule):
    "obtain cDNA sequences"
    def __init__(self, tm, srcDb, cdnaType):
        self.tm = tm
        ids = tm.getSrcSeqId(srcDb, cdnaType)
        fa = tm.getSrcFa(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType, "seq", accFile=FileIn(ids))
        cmd2 = ("awk", "-v", "srcDb="+srcDb.name, '''/^>/{$0 = ">" srcDb ":" substr($0, 2)} {print $0}''')
        CmdRule.__init__(self, Cmd((cmd1, cmd2), stdout=fa))

class CDnaMeta(CmdRule):
    "obtain meta data (CDS, geneName)"
    def __init__(self, tm, srcDb, cdnaType):
        "dbs are GenomeDef objects"
        ids = tm.getSrcSeqId(srcDb, cdnaType)
        meta = tm.getSrcMeta(srcDb, cdnaType)
        cmd1 = mkGbGetSeqsCmd(tm, srcDb, cdnaType, "ra")
        cmd2 = ("gbRaToMeta", srcDb.name)
        CmdRule.__init__(self, Cmd((cmd1, cmd2), stdout=meta))

class CDnaPolyA(CmdRule):
    "obtain poly information on cDNA sequences"
    # not currently used
    def __init__(self, tm, srcDb, cdnaType):
        self.tm = tm
        fa = tm.getSrcFa(srcDb, cdnaType)
        polya = tm.getSrcPoly(srcDb, cdnaType)
        cmd = ["faPolyASizes", FileIn(fa), FileOut(polya)]
        CmdRule.__init__(self, Cmd(cmd))

class UcscGenesAligns(CmdRule):
    "obtain PSLs and CDS files for UCSC genes"
    def __init__(self, tm, srcDb):
        self.tm = tm
        psl = tm.getSrcPsl(srcDb, CDnaType.ucscGenes)
        cds = tm.getSrcCds(srcDb, CDnaType.ucscGenes)

        cols = makeColumnWithSrcDbConcat(srcDb.name, "name") + ",chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds"
        cmd1 = ("hgsql", "-Ne",
                "select " + cols + " from knownGene",
                srcDb.name)
        cmd2 = ("genePredToFakePsl", srcDb.name, "/dev/stdin", "/dev/stdout", FileOut(cds))
        cmd3 = ("pslQueryUniq",)
        CmdRule.__init__(self, Cmd((cmd1, cmd2, cmd3), stdout=psl))
    
class UcscGenesMeta(CmdRule):
    "obtain meta files for UCSC genes, joining with CDS"
    def __init__(self, tm, srcDb):
        self.tm = tm
        cds = tm.getSrcCds(srcDb, CDnaType.ucscGenes)
        meta = tm.getSrcMeta(srcDb, CDnaType.ucscGenes)
        
        cmd1 = ("getUcscGenesMeta", srcDb.name, FileIn(cds), FileOut(meta))
        CmdRule.__init__(self, Cmd((cmd1,)))
    
class UcscGenesSeqs(CmdRule):
    def __init__(self, tm, srcDb):
        # warning: must get sequences from genePred, or it will not
        # match sizes in generated PSL.
        self.tm = tm
        fa = tm.getSrcFa(srcDb, CDnaType.ucscGenes)
        cmd1 = ("getGenePredFa", srcDb.name, "knownGene", FileOut(fa))
        CmdRule.__init__(self, Cmd((cmd1,)))

def getTablesLike(db, pat):
    return procOps.callProcLines(("hgsql", "-Ne", 'show tables like "'+pat+'";', db))

def getGencodeVersion(db):
    "return latest gencode version, parse from table names, or None"
    tbls = getTablesLike(db, "wgEncodeGencodeComp%")
    if len(tbls) == 0:
        return None
    # get most current version number
    versions = [int(re.search("^[^0-9]+([0-9]+)$", tbl).group(1)) for tbl in tbls]
    versions.sort(reverse=True)
    version = str(versions[0])
    # if mouse, then add back in M
    if tbls[0].startswith("wgEncodeGencodeCompVM"):
        version = "M" + version
    return version

def getGencodeTables(db):
    "return latest GENCODE tables (comp, attr) for latest, or None"
    version = getGencodeVersion(db)
    if version == None:
        return None
    else:
        return ("wgEncodeGencodeCompV"+str(version), "wgEncodeGencodeAttrsV"+version)

def getEnsemblTables(db):
    "return tables (ensGene, ensemblToGeneName) for latest, or None"
    if len(getTablesLike(db, "ensGene")) == 0:
        return None
    else:
        return ("ensGene", "ensemblToGeneName")

def getGencodeEnsemblTables(db):
    tbls = getGencodeTables(db)
    if tbls == None:
        tbls = getEnsemblTables(db)
    return tbls

# FIXME: also in commands, move to library
def isGenePredExt(db, gpTbl):
    # does it contain the exonFrames field?
    return procOps.callProc(["hgsql", "-Ne", "describe " + gpTbl, db]).find("exonFrames") >= 0

class EnsemblAligns(CmdRule):
    "obtain PSLs and CDS files for Ensembl or GENCODE genes"
    def __init__(self, tm, srcDb):
        self.tm = tm
        psl = tm.getSrcPsl(srcDb, CDnaType.ensembl)
        cds = tm.getSrcCds(srcDb, CDnaType.ensembl)
        
        tbls = getGencodeEnsemblTables(srcDb.name)

        cols =  makeColumnWithSrcDbConcat(srcDb.name, "name") + ', chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds'
        if isGenePredExt(srcDb.name, tbls[0]):
            cols += ', score, name2, cdsStartStat, cdsEndStat, exonFrames'
        cmd1 = ("hgsql", "-Ne", "select " + cols + " from " + tbls[0], srcDb.name)
        cmd2 = ("genePredToFakePsl", srcDb.name, "/dev/stdin", "/dev/stdout", FileOut(cds))
        cmd3 = ("pslQueryUniq",)
        CmdRule.__init__(self, Cmd((cmd1, cmd2, cmd3), stdout=psl))
    
class EnsemblMeta(CmdRule):
    "obtain meta files for Ensembl or GENCODE genes, joining with CDS"
    def __init__(self, tm, srcDb):
        self.tm = tm
        cds = tm.getSrcCds(srcDb, CDnaType.ensembl)
        meta = tm.getSrcMeta(srcDb, CDnaType.ensembl)
        
        gencodeTbls = getGencodeTables(srcDb.name)
        if gencodeTbls != None:
            cmd1 = ("getGencodeMeta", srcDb.name, gencodeTbls[1], FileIn(cds), FileOut(meta))
        else:
            cmd1 = ("getEnsemblMeta", srcDb.name, FileIn(cds), FileOut(meta))
        CmdRule.__init__(self, Cmd((cmd1,)))
    
class EnsemblSeqs(CmdRule):
    def __init__(self, tm, srcDb):
        # warning: must get sequences from genePred, or it will not
        # match sizes in generated PSL.
        self.tm = tm
        fa = tm.getSrcFa(srcDb, CDnaType.ensembl)
        
        geneTbl, ignore = getGencodeEnsemblTables(srcDb.name)
        cmd1 = ("getGenePredFa", srcDb.name, geneTbl, FileOut(fa))
        CmdRule.__init__(self, Cmd((cmd1,)))

def __createCdnaTypeRules(tm, srcDb, cdnaType):
    if cdnaType == CDnaType.ucscGenes:
        tm.exrun.addRule(UcscGenesAligns(tm, srcDb))
        tm.exrun.addRule(UcscGenesMeta(tm, srcDb))
        tm.exrun.addRule(UcscGenesSeqs(tm, srcDb))
    elif cdnaType == CDnaType.ensembl:
        tm.exrun.addRule(EnsemblAligns(tm, srcDb))
        tm.exrun.addRule(EnsemblMeta(tm, srcDb))
        tm.exrun.addRule(EnsemblSeqs(tm, srcDb))
    else:
        tm.exrun.addRule(CDnaAligns(tm, srcDb, cdnaType))
        tm.exrun.addRule(CDnaSeqs(tm, srcDb, cdnaType))

    tm.exrun.addRule(CDnaAlignStats(tm, srcDb, cdnaType))
    tm.exrun.addRule(CDnaSeqIds(tm, srcDb, cdnaType))
    if cdnaType in (CDnaType.refSeq, CDnaType.mrna):
        tm.exrun.addRule(CDnaMeta(tm, srcDb, cdnaType))

def createRules(tm, srcDb, cdnaTypes):
    "create rules to get data files for srcDb"
    for cdnaType in srcDb.cdnaTypes & cdnaTypes:
        __createCdnaTypeRules(tm, srcDb, cdnaType)
