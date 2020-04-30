
from Bio.Nexus import Trees
from pycbio.sys import fileOps
from transMap import hgDbNameParse
from collections import namedtuple

# FIXME: not currently used, but might make a nice, general module of
# pycbio


class PhyloDistances(object):
    """Phylogenetic distances from multiple tree files.  This uses the
    pair from the newest databases and keeps distance only on the organism,
    not the specific database.  This is done with the hope that newer distances are
    more accurate."""

    # stores parsed database with parent for picking the newest distance
    __hgDbPhyloDistance = namedtuple("hgDbPhyloDistance", ("hgOrg1", "hgVer1", "hgOrg2", "hgVer2", "dist"))

    def __init__(self, treeFiles):
        # indexed by (db1, db2), with the names sorted
        self._distances = {}
        self._buildDistanceTable(treeFiles)

    def get(self, hgDb1, hgDb2):
        "return distance or None is not defined"
        hgOrg1, hgVer1 = hgDbNameParse(hgDb1)
        hgOrg2, hgVer2 = hgDbNameParse(hgDb2)
        if hgOrg1 < hgOrg2:
            entry = self._distances.get((hgOrg1, hgOrg2))
        else:
            entry = self._distances.get((hgOrg2, hgOrg1))
        return entry.dist if entry is not None else None

    def _buildDistanceTable(self, treeFiles):
        # indexed by (hgOrg1, hgOrg2) which are sorted
        for treeFile in treeFiles:
            self._processTree(PhyloDistances.__loadTree(treeFile))

    def _processTree(self, tree):
        leafs = [tree.node(i) for i in tree.get_terminals()]
        for l1 in leafs:
            for l2 in leafs:
                if l1 != l2:
                    self._processLeaves(l1.data.taxon, l2.data.taxon, tree.distance(l1.id, l2.id))

    def _processLeaves(self, hgDb1, hgDb2, dist):
        hgOrg1, hgVer1 = hgDbNameParse(hgDb1)
        hgOrg2, hgVer2 = hgDbNameParse(hgDb2)
        if hgOrg1 < hgOrg2:
            self._processDbPair(hgOrg1, hgVer1, hgOrg2, hgVer2, dist)
        else:
            self._processDbPair(hgOrg2, hgVer2, hgOrg1, hgVer1, dist)

    def _processDbPair(self, hgOrg1, hgVer1, hgOrg2, hgVer2, dist):
        entry = self._distances.get((hgOrg1, hgOrg2))
        if (entry is None) or self._isNewer(hgVer1, hgVer2, entry):
            self._distances[(hgOrg1, hgOrg2)] = self._hgDbPhyloDistance(hgOrg1, hgVer1, hgOrg2, hgVer2, dist)

    @staticmethod
    def _isNewer(hgVer1, hgVer2, entry):
        return (hgVer1 > entry.hgVer1) or (hgVer2 > entry.hgVer2)

    @staticmethod
    def _parsePhastConsMod(treeFile):
        strtree = None
        for line in fileOps.iterLines(treeFile):
            if line.startswith("TREE:"):
                strtree = line[6:]
        if strtree is None:
            raise Exception("TREE: not found in phastCons model: " + treeFile)
        return strtree

    @staticmethod
    def _loadTree(treeFile):
        """Load a tree file, can either be a phastCons model containing Newick
        (ends in .mod) or a file containing a Newick tree.
        """
        if treeFile.endswith(".mod"):
            strtree = PhyloDistances.__parsePhastConsMod(treeFile)
        else:
            strtree = "".join(fileOps.readFileLines(treeFile))
        return Trees.Tree(strtree)
