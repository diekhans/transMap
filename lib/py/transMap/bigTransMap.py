"""
support for build bigTransMap files
"""

def _mkCommaList(elems):
    return ",".join(str(e) for e in elems) + ","

def _pslToBedIdent(psl):
    return int(10 * psl.identity())

def bigTransMapMakeRec(srcDb, srcPsl, mappedPsl, sequence, metadata=None):
    """object used to build bigTransMap rows"""
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
           mappedPsl.qStrand(),  # strand
           mappedPsl.tEnd,  # thickStart
           mappedPsl.tEnd,  # thickEnd
           "",  # reserved
           len(mappedPsl.exons),  # blockCount
           _mkCommaList([len(e) for e in mappedPsl.exons]),  # blockSizes
           _mkCommaList([e.tStart for e in mappedPsl.exons]),  # chromStarts
           mappedPsl.qStart,  # oChromStart
           mappedPsl.qEnd,  # oChromEnd
           mappedPsl.getTStrand()  # oStrand
           mappedPsl.qSize, # oChromSize
           _mkCommaList([e.qStart for e in mappedPsl.exons]), # oChromStarts
           sequence.seq,  # oSequence
           cds,  # oCDS
           mappedPsl.tSize  # chromSize
           mappedPsl.match  # match
           mappedPsl.misMatch  # misMatch
           mappedPsl.repMatch  # repMatch
           mappedPsl.nCount  # nCount
           1,  # seqType 1=nucleotide
           srcDb,  # srcDb
           srcPsl.tName,  # srcChrom
           srcPsl.tStart,  # srcChromStart
           srcPsl.tEnd,  # srcChromEnd
           _pslToBedIdent(srcPsl),  # srcScore
           geneName, # geneName
           geneId # geneId
    ]
    return row
