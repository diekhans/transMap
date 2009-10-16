#include "common.h"
#include "lib/alignPair.h"
#include "psl.h"
#include "linefile.h"
#include "lib/mapInfo.h"

/* constructor */
struct alignPair *alignPairNew(struct psl *psl, struct mapInfo *mi) {
    struct alignPair *ap;
    AllocVar(ap);
    ap->psl = psl;
    ap->mi = mi;
    return ap;
}

/* destructor */
void alignPairFree(struct alignPair **apPtr) {
    struct alignPair *ap = *apPtr;
    if (ap != NULL) {
        pslFree(&ap->psl);
        mapInfoFree(&ap->mi);
        freeMem(ap->key);
        freeMem(ap);
        *apPtr = NULL;
    }
}

/* construct a key for an alignment pair */
char* alignPairMkKey(char *qName, int qStart, int qEnd,
                     char *tName, int tStart, int tEnd) {
    static char key[256];
    safef(key, sizeof(key), "%s\t%d\t%d\t%s\t%d\t%d",
          qName, qStart, qEnd, tName, tStart, tEnd);
    return cloneString(key);
}

/* construct a key for an alignment pair from a psl */
char* alignPairMkPslKey(struct psl *psl) {
    return alignPairMkKey(psl->qName, psl->qStart, psl->qEnd,
                          psl->tName, psl->tStart, psl->tEnd);
}

/* construct a key for an alignment pair from a mapInfo object */
char* alignPairMkMapInfoKey(struct mapInfo *mi) {
    return alignPairMkKey(mi->mappedQName, mi->mappedQStart, mi->mappedQEnd,
                          mi->mappedTName, mi->mappedTStart, mi->mappedTEnd);
}

/* compare key fields in two alignPair objects, also sorting my query order */
boolean alignPairKeysCmp(struct alignPair *a, struct alignPair *b) {
    int d = strcmp(a->psl->qName, b->psl->qName);
    if (d == 0) {
        d = a->psl->qStart - b->psl->qStart;
        if (d == 0) {
            d = a->psl->qEnd - b->psl->qEnd;
            if (d == 0) {
                d = strcmp(a->psl->tName, b->psl->tName);
                if (d == 0) {
                    d = a->psl->tStart - b->psl->tStart;
                    if (d == 0) {
                        d = a->psl->tEnd - b->psl->tEnd;
                    }
                }
            }
        }
    }
    return d;
}

/* compare key fields in psl and mapInfo objects */
static boolean keyFieldsEq(struct psl *psl, struct mapInfo *mi) {
    return sameString(psl->qName, mi->mappedQName)
        && (psl->qStart == mi->mappedQStart)
        && (psl->qEnd == mi->mappedQEnd)
        && sameString(psl->tName, mi->mappedTName)
        && (psl->tStart == mi->mappedTStart)
        && (psl->tEnd == mi->mappedTEnd);
}

/* load an alignPair object from two parallel lineFiles */
struct alignPair *alignPairLoad(struct lineFile *pslLf, struct lineFile *miLf) {
    struct psl *psl = pslNext(pslLf);
    struct mapInfo *mi = mapInfoNext(miLf);
    if ((psl == NULL) && (mi == NULL)) {
        return NULL;
    } else if ((psl == NULL) || (mi == NULL)) {
        errAbort("mismatched number of records in psl and mapInfo files %s and %s", pslLf->fileName, miLf->fileName);
        return NULL;
    } else if (!keyFieldsEq(psl, mi)) {
        errAbort("mismatched key fields in psl and mapInfo files %s:%d and %s:%d", pslLf->fileName, pslLf->lineIx, miLf->fileName, miLf->lineIx);
        return NULL;
    } else {
        return alignPairNew(psl, mi);
    }
}

/* load all alignPair objects from two parallel files */
struct alignPair *alignPairLoadAll(char *pslFile, char *miFile) {
    struct alignPair *aps = NULL, *ap;
    struct lineFile *pslLf = lineFileOpen(pslFile, TRUE);
    struct lineFile *miLf = lineFileOpen(miFile, TRUE);
    while ((ap = alignPairLoad(pslLf, miLf)) != NULL) {
        slAddHead(&aps, ap);
    }
    lineFileClose(&pslLf);
    lineFileClose(&miLf);
    slReverse(&aps);
    return aps;
}

/* write an alignPair object to two parallel files */
void alignPairWrite(struct alignPair *ap, FILE *pslFh, FILE *miFh) {
    pslTabOut(ap->psl, pslFh);
    mapInfoTabOut(ap->mi, miFh);
}

/* write all alignPair objects to two parallel files */
void alignPairWriteAll(struct alignPair *aps, char *pslFile, char *miFile) {
    FILE *pslFh = mustOpen(pslFile, "w");
    FILE *miFh = mustOpen(miFile, "w");
    fputs(mapInfoHdrs, miFh);
    fputc('\n', miFh);
    for (struct alignPair *ap = aps; ap != NULL; ap = ap->next) {
        alignPairWrite(ap, pslFh, miFh);
    }
    carefulClose(&pslFh);
    carefulClose(&miFh);
}

/* sort align pair objects, first by query, then by target; this is also key order */
int alignPairQuerySortCmp(const void *va, const void *vb) {
    return alignPairKeysCmp(*((struct alignPair **)va), *((struct alignPair **)vb));
}
