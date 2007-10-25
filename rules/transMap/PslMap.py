"""ExRun rule class to create to map alignments from one genome to another using
pslMap
"""

from pycbio.exrun.ExRun import CmdRule
from transMap.GenomeDefs import ChainSelect

class PslMap(CmdRule):
    "Map a PSL from one genome to another"

    def __init__(self, exRun, id, srcPsl, destPre, chainFile, netFile=None, chainType=ChainSelect.all,
                 swapChains=False):
        """nets are ignored if not needed."""
        if (chainType == ChainSelect.syn) and (nets == None):
            raise Exception("selecting synthenic chains requires nets")
        # input
        self.srcPsl = exRun.getFile(srcPsl)
        self.chains = exRun.getFile(chainFile)
        requires = [self.srcPsl, self.chains]
        self.nets = None
        if nets != None:
            self.nets = exRun.getFile(netFile)
            requires.append(self.nets)

        # output
        self.rawPsl = exRun.getFile(destPre + ".rawPsl.gz")
        self.rawMapInfo = exRun.getFile(destPre + ".rawMapinfo.gz")

        # options
        self.chainType = chainType
        self.swapChains = swapChains

        CmdRule.__init__(self, id, None,
                         produces=(self.rawPsl, self.rawMapinfo),
                         requires=requires)

    def _netFilteredCmd(self):
        "get pipeline for net filtering of chains"
        cmds = []
        cmd = ["netFilter"]
        if self.chainType == ChainSelect.syn:
            cmd.append("-syn")
        cmd.append(self.nets.path)
        cmds.append(cmd)

        # netChainSubset
        cmds.append(["netChainSubset", "-wholeChains", "stdin", self.chains.path, "stdout"])
        return cmds
                
    def run(self):
        """run rule"""
        if self.nets != None:
            cmds = self._netFilteredCmd()
            chFile = "stdin"
        else:
            cmds = []
            chFile = self.chains.path
        cmd = ["pslMap", "-chainMapFile"]
        if self.swapChains:
            cmd.append("-swapMap")
        cmd += ("-mapInfo="+self.rawMapinfo.getUncompressed(), self.srcPsl.path,
                chFile, "stdout")
        cmds.append(cmd)
        cmds.append(("bzip2", "-c"))
        self.addCmd(Cmd(cmds, stdout=self.rawPsl.getTmpOut()))
        self.addCmd(("bzip2", "-f", self.rawMapinfo.path))
        CmdRule.run(self)
        self.rawPsl.installTmpOut()
