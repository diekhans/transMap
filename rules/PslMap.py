"""ExRun rule class to create to map alignments from one genome to another using
pslMap
"""

from pycbio.exrun.ExRun import CmdRule, Cmd, IFileRef, OFileRef
from transMap.GenomeDefs import ChainType, ChainOrient

# N.B. /dev/stdout is used instead of stdout, as it appears to fasta corruption problem

class PslMap(CmdRule):
    "Map a PSL from one genome to another"

    def __init__(self, tm, chains, cdnaType):
        """chains is GenomeDefs.Chains object"""
        #synFilter = (chains.srcDb.clade == chains.destDb.clade)
        synFilter = True  # alway synfilter

        # input
        chainsFile = tm.exrun.getFile(chains.getChainFile())
        srcPsl = tm.getSrcPsl(chains.srcDb, cdnaType)

        # output
        rawPsl = tm.getMappedRawPsl(chains.destDb, chains.srcDb, cdnaType)
        rawMapinfo = tm.getMappedRawMapinfo(chains.destDb, chains.srcDb, cdnaType)

        # build cmd
        cmds = []
        if synFilter:
            netsFile = tm.exrun.getFile(chains.getNetFile())
            cmds.append(["netFilter", "-syn", IFileRef(netsFile)])
            cmds.append(["netChainSubset", "-wholeChains", "-verbose=0", "/dev/stdin", IFileRef(chainsFile), "/dev/stdout"])
        cmd = ["pslMap", "-chainMapFile"]
        if chains.orient == ChainOrient.rev:
            cmd.append("-swapMap")
        cmd.append("-mapInfo="+OFileRef(rawMapinfo))
        cmd.append(IFileRef(srcPsl))
        if synFilter:
            cmd.append("stdin")
        else:
            cmd.append(IFileRef(chainsFile))
        cmd.append(OFileRef(rawPsl))
        cmds.append(cmd)

        CmdRule.__init__(self, Cmd(cmds))

class PslCDnaFilter(CmdRule):
    "filter raw psls"

    def __init__(self, tm, chains, cdnaType):
        # not using polyA, doesn't seem worth the cost
        
        # finished targets are more strigent
        if chains.destDb.isFinished():
            nearBest = "-globalNearBest=0.005"
        else:
            nearBest = "-globalNearBest=0.01"

        rawPsl = tm.getMappedRawPsl(chains.destDb, chains.srcDb, cdnaType)
        psl = tm.getMappedPsl(chains.destDb, chains.srcDb, cdnaType)
        filtStats = tm.getMappedFiltStats(chains.destDb, chains.srcDb, cdnaType)

        cmd1 = ("pslCDnaFilter", "-minQSize=20", "-minCover=0.20", "-bestOverlap",
                nearBest, IFileRef(rawPsl), OFileRef(psl))
        cmd2 = ("bin/pslQueryUniq",)
        CmdRule.__init__(self, Cmd((cmd1, cmd2), stderr=OFileRef(filtStats)))

class PslStats(CmdRule):
    "run pslStats"
    def __init__(self, tm, chains, cdnaType):
        psl = tm.getMappedPsl(chains.destDb, chains.srcDb, cdnaType)
        pslStats = tm.getMappedPslStats(chains.destDb, chains.srcDb, cdnaType)
        cmd = ("pslStats", IFileRef(psl), OFileRef(pslStats))
        CmdRule.__init__(self, Cmd(cmd))

def createRules(tm, chains):
    """create rules to generate mapping using chains for all cDNA types
    in the srcDb"""
    for cdnaType in chains.srcDb.cdnaTypes:
        tm.exrun.addRule(PslMap(tm, chains, cdnaType))
        tm.exrun.addRule(PslCDnaFilter(tm, chains, cdnaType))
        tm.exrun.addRule(PslStats(tm, chains, cdnaType))
