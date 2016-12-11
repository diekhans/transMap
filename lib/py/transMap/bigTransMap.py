"""
support for build bigTransMap files
"""
import os
from pycbio.sys import fileOps
from transMap import getSortProg
import pipettor


def _mkCommaList(elems):
    return ",".join(str(e) for e in elems) + ","


def _mkRelCommaList(off, elems):
    return ",".join(str(e - off) for e in elems) + ","


def _pslToBedIdent(psl):
    return int(10 * psl.identity())


def bigTransMapMakeRec(srcDb, srcPsl, mappedPsl, sequence, chainType, metadata):
    """object used to build bigTransMap rows; metadata maybe None"""
    if mappedPsl.getTStrand() == '-':
        mappedPsl = mappedPsl.reverseComplement()
    cds = geneName = geneId = ""
    if metadata is not None:
        if metadata.cds is not None:
            cds = metadata.cds
        if metadata.geneName is not None:
            geneName = metadata.geneName
        if metadata.geneId is not None:
            geneName = metadata.geneId

    row = [mappedPsl.tName,  # chrom
           mappedPsl.tStart,  # chromStart
           mappedPsl.tEnd,  # chromEnd
           mappedPsl.qName,  # name
           _pslToBedIdent(mappedPsl),  # score
           mappedPsl.getQStrand(),  # strand
           mappedPsl.tEnd,  # thickStart
           mappedPsl.tEnd,  # thickEnd
           0,  # reserved
           len(mappedPsl.blocks),  # blockCount
           _mkCommaList([len(e) for e in mappedPsl.blocks]),  # blockSizes
           _mkRelCommaList(mappedPsl.tStart, [e.tStart for e in mappedPsl.blocks]),  # chromStarts
           mappedPsl.qStart,  # oChromStart
           mappedPsl.qEnd,  # oChromEnd
           mappedPsl.getTStrand(),  # oStrand
           mappedPsl.qSize,  # oChromSize
           _mkCommaList([e.qStart for e in mappedPsl.blocks]),  # oChromStarts
           sequence.seq,  # oSequence
           cds,  # oCDS
           mappedPsl.tSize,  # chromSize
           mappedPsl.match,  # match
           mappedPsl.misMatch,  # misMatch
           mappedPsl.repMatch,  # repMatch
           mappedPsl.nCount,  # nCount
           1,  # seqType 1=nucleotide
           srcDb,  # srcDb
           srcPsl.tName,  # srcChrom
           srcPsl.tStart,  # srcChromStart
           srcPsl.tEnd,  # srcChromEnd
           _pslToBedIdent(srcPsl),  # srcScore
           geneName,
           geneId,
           chainType]
    return row


def _writeSortList(inFiles):
    tmpSortList = fileOps.tmpFileGet("transMap.merge-sort.")
    with open(tmpSortList, "w") as fh:
        for inFile in inFiles:
            fh.write(inFile)
            fh.write("\0")
    return tmpSortList


def _bigTransMapSortMerge(inFiles, outPreBigPsl):
    tmpSortList = _writeSortList(inFiles)
    outPreBigPslTmp = fileOps.atomicTmpFile(outPreBigPsl)
    sortCmd = [getSortProg(), "-k1,1", "-k2,2n", "-k3,3n", "-k4,4",
               "--merge", "--files0-from={}".format(tmpSortList),
               "--output={}".format(outPreBigPslTmp)]
    pipettor.run(sortCmd)
    fileOps.atomicInstall(outPreBigPslTmp, outPreBigPsl)
    os.unlink(tmpSortList)


def bigTransMapSortMerge(inFiles, outPreBigPsl):
    # sort generates error on empty input, so handle this case
    fileOps.ensureFileDir(outPreBigPsl)
    if len(inFiles) == 0:
        open(outPreBigPsl, "w").close()  # create empty
    else:
        _bigTransMapSortMerge(inFiles, outPreBigPsl)
