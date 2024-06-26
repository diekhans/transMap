"""
support for build bigTransMap files
"""
from builtins import str
import os
from pycbio.sys import fileOps
from transMap import getSortProg, alignIdToSrcId, srcIdToAccv
from collections import namedtuple
import pipettor


def _toIntArray(a):
    "ensure values in array are ints"
    return [int(v) for v in a]


class BigTransMap(namedtuple("BigTransMap",
                             ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
                              "thickEnd", "reserved", "blockCount", "blockSizes", "chromStarts",
                              "oChromStart", "oChromEnd", "oStrand", "oChromSize", "oChromStarts",
                              "oSequence", "oCDS", "chromSize", "match", "misMatch", "repMatch", "nCount",
                              "seqType", "srcDb", "srcTransId", "srcChrom", "srcChromStart", "srcChromEnd",
                              "srcIdent", "srcAligned", "geneName", "geneId", "geneType", "transcriptType",
                              "chainType", "commonName", "scientificName", "orgAbbrev"))):
    __slots__ = ()

    def __init__(self, chrom, chromStart, chromEnd, name, score, strand,
                 thickStart, thickEnd, reserved, blockCount, blockSizes,
                 chromStarts, oChromStart, oChromEnd, oStrand, oChromSize,
                 oChromStarts, oSequence, oCDS, chromSize, match, misMatch,
                 repMatch, nCount, seqType, srcDb, srcTransId, srcChrom,
                 srcChromStart, srcChromEnd, srcIdent, srcAligned, geneName,
                 geneId, geneType, transcriptType, chainType, commonName,
                 scientificName, orgAbbrev):
        super(BigTransMap, self).__init__(chrom, int(chromStart),
                                          int(chromEnd), name, int(score),
                                          strand, int(thickStart),
                                          int(thickEnd), int(reserved),
                                          int(blockCount),
                                          _toIntArray(blockSizes),
                                          _toIntArray(chromStarts),
                                          int(oChromStart), int(oChromEnd),
                                          oStrand, int(oChromSize),
                                          _toIntArray(oChromStarts),
                                          oSequence, oCDS, int(chromSize),
                                          int(match), int(misMatch),
                                          int(repMatch), int(nCount),
                                          int(seqType), srcDb, srcTransId,
                                          srcChrom, int(srcChromStart),
                                          int(srcChromEnd), int(srcIdent),
                                          int(srcAligned), geneName, geneId,
                                          geneType, transcriptType, chainType,
                                          commonName, scientificName,
                                          orgAbbrev)


def _mkCommaList(elems):
    return ",".join(str(e) for e in elems) + ","


def _mkRelCommaList(off, elems):
    return ",".join(str(e - off) for e in elems) + ","


def _pslToBedIdent(psl):
    return int(1000 * psl.identity())


def _pslToBedAligned(psl):
    return int(1000 * psl.queryAligned())


def bigTransMapMakeRec(srcDb, srcPsl, mappedPsl, sequence, chainType, metadata,
                       commonName, scientificName, orgAbbrev):
    """object used to build bigTransMap rows; metadata maybe None"""
    if mappedPsl.getTStrand() == '-':
        mappedPsl = mappedPsl.reverseComplement()
    cds = geneName = geneId = geneType = transcriptType = ""
    if metadata is not None:
        if metadata.cds is not None:
            cds = metadata.cds
        if metadata.geneName is not None:
            geneName = metadata.geneName
        if metadata.geneId is not None:
            geneId = metadata.geneId
        if metadata.geneType is not None:
            geneType = metadata.geneType
        if metadata.transcriptType is not None:
            transcriptType = metadata.transcriptType

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
           srcIdToAccv(alignIdToSrcId(srcPsl.qName)),  # srcTransId
           srcPsl.tName,  # srcChrom
           srcPsl.tStart,  # srcChromStart
           srcPsl.tEnd,  # srcChromEnd
           _pslToBedIdent(srcPsl),  # srcIdent
           _pslToBedAligned(srcPsl),  # srcAligned
           geneName,
           geneId,
           geneType,
           transcriptType,
           chainType,
           commonName,
           scientificName,
           orgAbbrev]
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
    with fileOps.AtomicFileCreate(outPreBigPsl) as outPreBigPslTmp:
        sortCmd = [getSortProg(), "-k1,1", "-k2,2n", "-k3,3n", "-k4,4",
                   "--merge", "--files0-from={}".format(tmpSortList),
                   "--output={}".format(outPreBigPslTmp)]
        pipettor.run(sortCmd)
    os.unlink(tmpSortList)


def bigTransMapSortMerge(inFiles, outPreBigPsl):
    # sort generates error on empty input, so handle this case
    fileOps.ensureFileDir(outPreBigPsl)
    if len(inFiles) == 0:
        open(outPreBigPsl, "w").close()  # create empty
    else:
        _bigTransMapSortMerge(inFiles, outPreBigPsl)
