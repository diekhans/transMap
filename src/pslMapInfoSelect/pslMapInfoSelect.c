static const char *usageMsg =
    "pslMapInfoSelect mapInfoIn selectPsl mapInfoOut\n"
    "\n"
    "Select the subset of mapInfo records that are in selectPsl.\n"
    "This also updates the mappedQName record to match the uniqueness\n"
    "suffix that was add to qName.  The mapInfoIn chainSub column is\n"
    "optional.\n"
    "\n"
    "Options:\n"
    "   -dropUniq2nd - drop 2nd uniqueness suffix\n";

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "localmem.h"
#include "psl.h"
#include "options.h"
#include "lib/mapInfo.h"
#include "lib/alignPair.h"

/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"dropUniq2nd", OPTION_BOOLEAN},
    {NULL, 0}
};
static boolean gDropUniq2nd;

/* entry in the select tbl */
struct selectEntry {
    char *mappedQName;
    boolean done;
};

/* drop second uniqueness suffix in a psl qName.  
 * NM_205019.1-1.1 -> NM_205019.1-1
 * WARNING: static return. */
static char *dropUniq2nd(char *qName) {
    static char buf[64];
    safecpy(buf, sizeof(buf), qName);
    char *dash = strrchr(buf, '-');
    char *dot = strrchr(buf, '.');
    if ((dash == NULL) || (dot == NULL) || (dot < dash)) {
        errAbort("can't parse uniqueness suffix in %s", qName);
    }
    *dot = '\0';
    return buf;
}

/* add a psl to the select table. */
static void addPsl(struct hash *selectTbl,
                   struct psl *psl) {
    char *srcQName = psl->qName;
    if (gDropUniq2nd) {
        srcQName = dropUniq2nd(srcQName);
    }
    struct selectEntry *se;
    lmAllocVar(selectTbl->lm, se);
    se->mappedQName = lmCloneString(selectTbl->lm, psl->qName);
    char* key = alignPairMkKey(srcQName, psl->qStart, psl->qEnd,
                               psl->tName, psl->tStart, psl->tEnd);
    hashAdd(selectTbl, key, se);
    freeMem(key);
}

/* generate select table from psls */
static struct hash* buildSelectTbl(char *selectPsl) {
    struct hash *selectTbl = hashNew(24);
    struct lineFile *lf = pslFileOpen(selectPsl);
    struct psl *psl;

    while ((psl = pslNext(lf)) != NULL) {
        addPsl(selectTbl, psl);
        pslFree(&psl);
    }
    lineFileClose(&lf);
    return selectTbl;
}

/* check if mapInfo row be kept, returning selectEntry or NULL */
static struct selectEntry *findSelectEntry(struct hash *selectTbl,
                                           struct mapInfo *mi) {
    char *key = alignPairMkMapInfoKey(mi);
    struct selectEntry *se =hashFindVal(selectTbl, key);
    if ((se != NULL) && se->done) {
        errAbort("select entry already used: %s", key); 
    }
    freeMem(key);
    return se;
}

/* process one mapInfo row */
static void processMapInfo(struct hash *selectTbl,
                           struct lineFile *inLf,
                           FILE *outFh,
                           struct mapInfo *mi) {
    struct selectEntry *se = findSelectEntry(selectTbl, mi);
    if (se != NULL) {
        char *hold = mi->mappedQName;
        mi->mappedQName = se->mappedQName;
        mapInfoTabOut(mi, outFh);
        mi->mappedQName = hold;
        se->done = TRUE;
    }
}
                          

/* filter mapInfo using selected PSLs */
static void pslMapInfoSelect(char *mapInfoIn,
                             char *selectPsl,
                             char *mapInfoOut) {
    struct hash *selectTbl = buildSelectTbl(selectPsl);
    struct lineFile *inLf = lineFileOpen(mapInfoIn, TRUE);
    FILE *outFh = mustOpen(mapInfoOut, "w");
    fputs(mapInfoHdrs, outFh);
    fputc('\n', outFh);

    struct mapInfo *mi;
    while ((mi = mapInfoNext(inLf)) != NULL) {
        processMapInfo(selectTbl, inLf, outFh, mi);
        mapInfoFree(&mi);
    }
    carefulClose(&outFh);
    lineFileClose(&inLf);
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    gDropUniq2nd = optionExists("dropUniq2nd");
    pslMapInfoSelect(argv[1], argv[2], argv[3]);
    return 0;
}

