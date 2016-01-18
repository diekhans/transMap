
from Bio.Nexus import Trees
from pycbio.sys import fileOps


def _parsePhastConsMod(treeFile):
    strtree = None
    for line in fileOps.iterLines(treeFile):
        if line.startswith("TREE:"):
            strtree = line[6:]
    if strtree == None:
        raise Exception("TREE: not found in phastCons model: " + treeFile)
    return strtree

def load(treeFile):
    """Load a tree file, can either be a phastCons model containing Newick
    (ends in .mod) or a file containing a Newick tree.
    """
    if treeFile.endswith(".mod"):
        strtree = _parsePhastConsMod(treeFile)
    else:
        strtree = "".join(fileOps.readFileLines(treeFile))
    return Trees.Tree(strtree)

