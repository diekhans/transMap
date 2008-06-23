static const char *usageMsg =
    "pslFaPartition targetSize inPsl inFa outRootDir\n"
    "\n"
    "Partition psl and fa files in whole sequences so that they have\n"
    "approximately targetSize bases per output file.  Files will be parallel,\n"
    "each containing the same sequences.\n"
    "\n"
    "Directory must be empty, as files maybe appended\n"
    "\n"
    "Options:\n"
    "  -prefixList=out - write list of prefixes (less .psl or .fa) create to this file.\n";

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "sqlNum.h"
#include "psl.h"
#include "fa.h"
#include "hash.h"
#include "portable.h"
#include "localmem.h"


/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"prefixList", OPTION_STRING},
    {NULL, 0},
};

/* partition table */
struct partTbl {
    struct hash *idToPart;
    unsigned nextPartNum;
    char *curPart;          // string number for current partition
    FILE *partFh;           // currently open file for curPart.
};

/* create partition table */
static struct partTbl *partTblNew() {
    struct partTbl *pt;
    AllocVar(pt);
    pt->nextPartNum = 0;
    pt->idToPart = hashNew(25);
    return pt;
}

/* move to the next partition */
static void partTblNext(struct partTbl *pt) {
    char curPart[64];
    safef(curPart, sizeof(curPart), "%05d", pt->nextPartNum++);
    pt->curPart = lmCloneString(pt->idToPart->lm, curPart);
}

/* Add an entry to the partition table.  This requires that all ids are
 * unique. */
static void partTblAdd(struct partTbl *pt, char *id) {
    struct hashEl *hel = hashStore(pt->idToPart, id);
    if (hel != NULL) {
        errAbort("duplicate id: %s", id);
    }
    hel->val = pt->curPart;
}

/* close currently open file. */
static void partTblFileClose(struct partTbl *pt) {
    carefulClose(&pt->partFh);
    pt->curPart = NULL;
}

/* Generate file path for a partition. WARNING static return */
static char *partTblFilePath(struct partTbl *pt, char *outDir, char *part, char *ext) {
    static char path[PATH_LEN];
    safef(path, sizeof(path), "%s/%s.%s", outDir, part, ext);
    return path;
}

/* Open the next partition file for writing */
static void partTblFileNext(struct partTbl *pt, char *outDir, char *ext) {
    partTblFileClose(pt);
    partTblNext(pt);
    pt->partFh = mustOpen(partTblFilePath(pt, outDir, pt->curPart, ext), "w");
}

/* Get a partition file for appending by id, caching last open.  Return NULL
 * if id is not found. */
static FILE *partTblFileGet(struct partTbl *pt, char *id, char *outDir, char *ext) {
    char *part = hashFindVal(pt->idToPart, id);
    if (part == NULL) {
        return NULL;
    }
    if (part != pt->curPart) {
        partTblFileClose(pt);
        pt->curPart = part;
        pt->partFh = mustOpen(partTblFilePath(pt, outDir, pt->curPart, ext), "a");
    }
    return pt->partFh;
}

/* Partition psl */
static struct partTbl *pslPartition(unsigned targetSize, char *inPslFile, char *outDir) {
    struct partTbl *pt = partTblNew();
    struct lineFile *inLf = pslFileOpen(inPslFile);
    unsigned curSize = 0;
    struct psl *psl;
    while ((psl = pslNext(inLf)) != NULL) {
        if (pt->partFh == NULL) {
            partTblFileNext(pt, outDir, "psl");
        }
        pslTabOut(psl, pt->partFh);
        partTblAdd(pt, psl->qName);
        curSize += psl->qSize;
        if (curSize >= targetSize) {
            partTblFileClose(pt);
            curSize = 0;
        }
        pslFree(&psl);
    }
    partTblFileClose(pt);
    lineFileClose(&inLf);
    return pt;
}

/* partition one sequence */
static void faPartitionSeq(struct partTbl *pt, char *outDir, DNA *seq, int size, char *header) {
    char buf[256];
    safecpy(buf, sizeof(buf), header);
    char *id = firstWordInLine(buf);
    FILE *fh = partTblFileGet(pt, id, outDir, "fa");
    if (fh != NULL) {
        faWriteNext(fh, header, seq, size);
    }
}

/* Partition fasta */
static void faPartition(struct partTbl *pt, char *inFaFile, char *outDir) {
    struct lineFile *inLf = lineFileOpen(inFaFile, TRUE);
    DNA *seq;
    int size;
    char *header;
    while (faSpeedReadNext(inLf, &seq, &size, &header)) {
        faPartitionSeq(pt, outDir, seq, size, header);
    }
    faFreeFastBuf();
    lineFileClose(&inLf);
}

/* write prefixes to a file */
static void listPrefixes(struct partTbl *pt, char *outDir, char *prefixList) {
    FILE *fh = mustOpen(prefixList, "w");
    int iPart;
    for (iPart = 0; iPart < pt->nextPartNum; iPart++) {
        fprintf(fh, "%s/%05d\n", outDir, iPart);
    }
    carefulClose(&fh);
}


/* partition psl and fa */
static void pslFaPartition(unsigned targetSize, char *inPslFile, char *inFaFile,
                           char *outDir, char *prefixList) {
    if (listDir(outDir, "*") != NULL) {
        errAbort("output directory not empty: %s", outDir);
    }
    struct partTbl *pt = pslPartition(targetSize, inPslFile, outDir);
    faPartition(pt, inFaFile, outDir);
    if (prefixList != NULL) {
        listPrefixes(pt, outDir, prefixList);
    }
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    pslFaPartition(sqlUnsigned(argv[1]), argv[1], argv[2], argv[3],
                   optionVal("prefixList", NULL));
    return 0;
}

