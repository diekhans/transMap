static const char *usageMsg =
    "pslFaPartition targetSize inPsl inFa outRootDir\n"
    "\n"
    "Partition psl and fa files in whole sequences so that they have\n"
    "approximately targetSize bases per output file.  Files will be parallel,\n"
    "each containing the same sequences.\n"
    "\n"
    "Options:\n"
    "  -partList=out - write list of part numbers created to this file.\n";

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
    {"partList", OPTION_STRING},
    {NULL, 0},
};

/* chop uniq suffix from acc. WARNING: static return */
static char *idToAccv(char *id) {
    static char accv[64];
    safecpy(accv, sizeof(accv), id);
    char *dash = strrchr(accv, '-');
    if (dash != NULL) {
        *dash = '\0';
    }
    return accv;
}

/* partition table */
struct partTbl {
    struct hash *accvToParts; // indexed by accv, value is slName of parts
    unsigned nextPartNum;
    char *curPart;            // string number for current partition
    FILE *partFh;             // currently open file for curPart.
};

/* create partition table */
static struct partTbl *partTblNew() {
    struct partTbl *pt;
    AllocVar(pt);
    pt->nextPartNum = 0;
    pt->accvToParts = hashNew(23);
    return pt;
}

/* move to the next partition */
static void partTblNext(struct partTbl *pt) {
    char curPart[64];
    safef(curPart, sizeof(curPart), "%05d", pt->nextPartNum++);
    pt->curPart = lmCloneString(pt->accvToParts->lm, curPart);
}

/* Add an entry to the partition table. don't duplicate if in existing
 * partition */
static void partTblAdd(struct partTbl *pt, char *id) {
    char *accv = idToAccv(id);
    struct hashEl *hel = hashStore(pt->accvToParts, accv);
    struct slName **parts = (struct slName**)&hel->val;
    if (!slNameInList(*parts, pt->curPart)) {
        slSafeAddHead(parts, lmSlName(pt->accvToParts->lm, pt->curPart));
    }
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

/* Get list of partitions for an accv, NULL if none */
static struct slName *partTblGetAccvParts(struct partTbl *pt, char *accv) {
    return hashFindVal(pt->accvToParts, accv);
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

/* table of open fasta files */
struct faPart {
    int numOpen;
    struct hash *partToFh;
};

/* build/rebuild hash */
static struct faPart *faPartNew() {
    struct faPart *faPart;
    AllocVar(faPart);
    faPart->partToFh = hashNew(20);
    return faPart;
}

/* close open fasta files */
static  void faPartClose(struct faPart *faPart) {
    struct hashCookie cookie = hashFirst(faPart->partToFh);
    struct hashEl *hel;
    while ((hel = hashNext(&cookie)) != NULL) {
        carefulClose((FILE **)&(hel->val));
    }
    faPart->numOpen = 0;
}

/* get fasta FILE for a partition */
static FILE *faPartGetFh(struct partTbl *pt, struct faPart *faPart, char *outDir, char *part) {
    struct hashEl *hel = hashStore(faPart->partToFh, part);
    if (hel->val == NULL) {
        if (faPart->numOpen >= 1000) {
            faPartClose(faPart);
        }
        hel->val = mustOpen(partTblFilePath(pt, outDir, part, "fa"), "a");
        faPart->numOpen++;
    }
    return hel->val;
}

/* partition one sequence */
static void faPartitionSeq(struct partTbl *pt, struct faPart *faPart, char *outDir, DNA *seq, int size, char *header) {
    char buf[256];
    safecpy(buf, sizeof(buf), header);
    char *accv = firstWordInLine(buf);
    struct slName *part, *parts = partTblGetAccvParts(pt, accv);
    for (part = parts; part != NULL; part = part->next) {
        FILE *fh = faPartGetFh(pt, faPart, outDir, part->name);
        faWriteNext(fh, header, seq, size);
    }
}

/* Partition fasta */
static void faPartition(struct partTbl *pt, char *inFaFile, char *outDir) {
    struct faPart *faPart = faPartNew();
    struct lineFile *inLf = lineFileOpen(inFaFile, TRUE);
    DNA *seq;
    int size;
    char *header;
    while (faSpeedReadNext(inLf, &seq, &size, &header)) {
        faPartitionSeq(pt, faPart, outDir, seq, size, header);
    }
    faFreeFastBuf();
    lineFileClose(&inLf);
    faPartClose(faPart);
}

/* write part prefixes to a file */
static void listParts(struct partTbl *pt, char *partList) {
    FILE *fh = mustOpen(partList, "w");
    int iPart;
    for (iPart = 0; iPart < pt->nextPartNum; iPart++) {
        fprintf(fh, "%05d\n", iPart);
    }
    carefulClose(&fh);
}

/* partition psl and fa */
static void pslFaPartition(unsigned targetSize, char *inPslFile, char *inFaFile,
                           char *outDir, char *partList) {
    if (listDir(outDir, "*") != NULL) {
        errAbort("output directory not empty: %s", outDir);
    }
    struct partTbl *pt = pslPartition(targetSize, inPslFile, outDir);
    faPartition(pt, inFaFile, outDir);
    if (partList != NULL) {
        listParts(pt, partList);
    }
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    pslFaPartition(sqlUnsigned(argv[1]), argv[2], argv[3], argv[4],
                   optionVal("partList", NULL));
    return 0;
}

