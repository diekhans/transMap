import os
import sys
import glob
from transMap.genomeData import ChainType, Chains

# FIXME: no data, this should not be a class

debug = False


class ChainsFinder(object):
    """Locate all chains between active species and find branch lengths"""

    # /hive/data/genomes/hg18/bed/blastz.calJac1/axtChain/
    #     hg18.calJac1.all.chain.gz
    #     hg18.calJac1.net.gz
    #     hg18.calJac1.rbest.chain.gz
    #     hg18.calJac1.rbest.net.gz
    # also lastz directories

    def findChains(self, hgDbs):
        """get genomeData.Chains for all pairings between provided databases"""
        chains = []

        # look for all combinations in both directions
        for hgDb1 in hgDbs:
            for hgDb2 in hgDbs:
                if hgDb1 != hgDb2:
                    self._findChains(hgDb1, hgDb2, chains)
        return chains

    def _findChains(self, destHgDb, srcHgDb, chains):
        chainsDirs = self._findChainsDirs(destHgDb, srcHgDb)
        self._addOneIfExists(destHgDb, srcHgDb, chainsDirs, ChainType.all, "all.chain", "net", chains)
        self._addOneIfExists(destHgDb, srcHgDb, chainsDirs, ChainType.rbest, "rbest.chain", "rbest.net", chains)
        # syntenic chains are generated from the all chains
        self._addOneIfExists(destHgDb, srcHgDb, chainsDirs, ChainType.syn, "all.chain", "net", chains)

    def _findChainsDirs(self, destHgDb, srcHgDb):
        "return list of possibles to chain directory, they might not actually contain chains"
        return self._findChainsDirsAligner(destHgDb, srcHgDb, "lastz") \
            + self._findChainsDirsAligner(destHgDb, srcHgDb, "blastz")

    def _findChainsDirsAligner(self, destHgDb, srcHgDb, aligner):
        """Return list of paths to possible chains dirs.  These might not actually
        contain the nets and chains.  It's simpler debug to get possible chain """
        chainsDirs = []
        base = self._chainsBasePath(destHgDb, srcHgDb, aligner)
        if os.path.exists(base):
            chainsDirs.append(base)
        chainsDir = base + ".swap"
        if os.path.exists(chainsDir):
            chainsDirs.append(chainsDir)
        # look for data extensions w/wo .swap
        matchDirs = glob.glob(base + ".*-*-*")
        matchDirs.sort(reverse=True)
        for chainsDir in matchDirs:
            if os.path.exists(chainsDir):
                chainsDirs.append(chainsDir)

        # check for missing symlink to names like blastzPanTro2.2006-03-28
        pat = self._chainsBasePath(destHgDb, srcHgDb[0].upper() + srcHgDb[1:], aligner) + ".*-*-*"
        matchDirs = glob.glob(pat)
        matchDirs.sort()
        for chainsDir in matchDirs:
            if os.path.exists(chainsDir):
                chainsDirs.append(chainsDir)
        return chainsDirs

    def _chainsBasePath(self, destHgDb, srcHgDb, meth):
        return os.path.join("/hive/data/genomes", destHgDb, "bed", "{}.{}".format(meth, srcHgDb))

    def _getLastzFile(self, destHgDb, srcHgDb, chainsDir, ext):
        """generate and check for one of the lastz alignment files, with and
        without compressed extensions.  Return path if one exists, otherwise
        None"""
        bzpre = os.path.join(chainsDir, "axtChain", "{}.{}.".format(destHgDb, srcHgDb))
        p = bzpre + ext + ".gz"  # common case
        if debug:
            sys.stderr.write("_getLastzFile: {} {} {}: {}\n".format(destHgDb, srcHgDb, os.path.exists(p), p))
        if os.path.exists(p):
            return p
        p = bzpre + ext
        if debug:
            sys.stderr.write("_getLastzFile: {} {} {}: {}\n".format(destHgDb, srcHgDb, os.path.exists(p), p))
        if os.path.exists(p):
            return p
        if debug:
            sys.stderr.write("_getLastzFile: {} {} {}\n".format(destHgDb, srcHgDb, "notFound"))
        return None

    def _addOneIfExists(self, destHgDb, srcHgDb, chainsDirs, chainType, chainExt, netExt, chains):
        "add first existing chain/net on list of chainsDirs"
        for chainsDir in chainsDirs:
            if self._addIfExists(destHgDb, srcHgDb, chainsDir, chainType, chainExt, netExt, chains):
                break

    def _addIfExists(self, destHgDb, srcHgDb, chainsDir, chainType, chainExt, netExt, chains):
        """add a set of nets/chains if they exist for the specified type"""
        chainFile = self._getLastzFile(destHgDb, srcHgDb, chainsDir, chainExt)
        if chainFile is not None:
            netFile = self._getLastzFile(destHgDb, srcHgDb, chainsDir, netExt)
            if netFile is not None:
                chains.append(Chains(srcHgDb, destHgDb, chainType, chainFile, netFile))
                return True
        return False
