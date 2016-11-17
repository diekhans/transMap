import os
import sys
import glob
from transMap.genomeData import ChainType, Chains


debug = False


class ChainsFinder(object):
    """Locate all chains between active species and find branch lengths"""

    # /hive/data/genomes/hg18/bed/blastz.calJac1/axtChain/
    #     hg18.calJac1.all.chain.gz
    #     hg18.calJac1.net.gz
    #     hg18.calJac1.rbest.chain.gz
    #     hg18.calJac1.rbest.net.gz
    # /hive/data/genomes/hg18/bed/blastz.ponAbe2/axtChain/
    #     hg18.ponAbe2.syn.net.gz
    # also lastz directories

    def __init__(self, distances):
        self.__distances = distances

    def findChains(self, hgDbs):
        """get genomeData.Chains for all pairings between provided databases"""
        chains = []

        # look for all combinations in both directions
        for hgDb1 in hgDbs:
            for hgDb2 in hgDbs:
                if hgDb1 != hgDb2:
                    self.__findChains(hgDb1, hgDb2, chains)
        return chains

    def __findChains(self, targetHgDb, queryHgDb, chains):
        chainsDirs = self.__findChainsDirs(targetHgDb, queryHgDb)
        self.__addOneIfExists(targetHgDb, queryHgDb, chainsDirs, ChainType.all, "all.chain", "net", chains)
        self.__addOneIfExists(targetHgDb, queryHgDb, chainsDirs, ChainType.syn, "all.chain", "syn.net", chains)
        self.__addOneIfExists(targetHgDb, queryHgDb, chainsDirs, ChainType.rbest, "rbest.chain", "rbest.net", chains)

    def __findChainsDirs(self, targetHgDb, queryHgDb):
        "return list of possibles to chain directory, they might not actually contain chains"
        return self.__findChainsDirsAligner(targetHgDb, queryHgDb, "lastz") \
            + self.__findChainsDirsAligner(targetHgDb, queryHgDb, "blastz")

    def __findChainsDirsAligner(self, targetHgDb, queryHgDb, aligner):
        """Return list of paths to possible chains dirs.  These might not actually
        contain the nets and chains.  It's simpler debug to get possible chain """
        chainsDirs = []
        base = self.__chainsBasePath(targetHgDb, queryHgDb, aligner)
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
        pat = self.__chainsBasePath(targetHgDb, queryHgDb[0].upper() + queryHgDb[1:], aligner) + ".*-*-*"
        matchDirs = glob.glob(pat)
        matchDirs.sort()
        for chainsDir in matchDirs:
            if os.path.exists(chainsDir):
                chainsDirs.append(chainsDir)
        return chainsDirs

    def __chainsBasePath(self, targetHgDb, queryHgDb, meth):
        return os.path.join("/hive/data/genomes", targetHgDb, "bed", "{}.{}".format(meth, queryHgDb))

    def __getLastzFile(self, targetHgDb, queryHgDb, chainsDir, ext):
        """generate and check for one of the blastz alignment files, with and
        without compressed extensions.  Return path if one exists, otherwise
        None"""
        bzpre = os.path.join(chainsDir, "axtChain", "{}.{}.".format(targetHgDb, queryHgDb))
        p = bzpre + ext + ".gz"  # common case
        if debug:
            sys.stderr.write("__getLastzFile: {} {} {}: {}\n".format(targetHgDb, queryHgDb, os.path.exists(p), p))
        if os.path.exists(p):
            return p
        p = bzpre + ext
        if debug:
            sys.stderr.write("__getLastzFile: {} {} {}: {}\n".format(targetHgDb, queryHgDb, os.path.exists(p), p))
        if os.path.exists(p):
            return p
        if debug:
            sys.stderr.write("__getLastzFile: {} {} {}\n".format(targetHgDb, queryHgDb, "notFound"))
        return None

    def __addOneIfExists(self, targetHgDb, queryHgDb, chainsDirs, chainType, chainExt, netExt, chains):
        "add first existing chain/net on list of chainsDirs"
        for chainsDir in chainsDirs:
            if self.__addIfExists(targetHgDb, queryHgDb, chainsDir, chainType, chainExt, netExt, chains):
                break

    def __addIfExists(self, targetHgDb, queryHgDb, chainsDir, chainType, chainExt, netExt, chains):
        """add a set of nets/chains if they exist for the specified type"""
        chainFile = self.__getLastzFile(targetHgDb, queryHgDb, chainsDir, chainExt)
        if chainFile is not None:
            netFile = self.__getLastzFile(targetHgDb, queryHgDb, chainsDir, netExt)
            if netFile is not None:
                chains.append(Chains(targetHgDb, queryHgDb, chainType, chainFile, netFile,
                                     self.__distances.get(targetHgDb, queryHgDb)))
                return True
        return False
