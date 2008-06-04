static const char *usageMsg =
    "augustusGenePredMk -cdsFile=cdsFile srcPsl mappedPsl augustusGp\n"
    "\n"
    "Generate genePreds for augustus which include a field indicating\n"
    "which gaps correspond to exons in the source\n";

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "psl.h"
#include "genePred.h"
#include "genbank.h"
#include "options.h"
#include "transMapInfo.h"
#include "lib/mapInfo.h"

/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"cdsFile", OPTION_STRING},
    {NULL, 0}
};

/* get the accession from a src or mapped id.  WARNING: static return */
static char *idToAcc(char *id) {
    static char buf[64];
    safecpy(buf, sizeof(buf), id);
    char *dash = strrchr(buf, '-');
    if (dash == NULL) {
        errAbort("invalid src or mapped id: %s", id);
    }
    *dash = '\0';
    return buf;
}

/* get the src from a mapped id.  WARNING: static return */
static char *mappedIdToSrcId(char *id) {
    static char buf[64];
    safecpy(buf, sizeof(buf), id);
    char *dash = strrchr(buf, '-');
    char *dot = strrchr(buf, '.');
    if ((dash == NULL) || (dot == NULL) || (dash > dot)) {
        errAbort("invalid mapped id: %s", id);
    }
    return buf;
}


/* load a hash by name of CDS */
static struct hash *cdsTblLoad(char *cdsFile) {
    struct hash *cdsTbl = hashNew(32*1048576);
    struct lineFile *inLf = lineFileOpen(cdsFile, TRUE);
    char *row[2];
    while (lineFileNextRowTab(inLf, row, ArraySize(row))) {
        struct genbankCds cds;
        if (genbankCdsParts(row[1], &cds)) {
            hashAdd(cdsTbl, row[0], lmCloneMem(cdsTbl->ln, cds, sizeof(cds)));
        }
    }
    lineFileClose(&inLf);
    return cdsTbl;
}

/* look up CDS, NULL if no CDS or table */
struct genbankCds *getCds(struct hash *cdsTbl,
                          char *id)

/* load a hash by name of source psls */
static struct hash *srcGeneLoad(struct hash *cdsTbl,
                                char *srcPslFile) {
    struct hash *srcGenes = hashNew(64*1048576);
    struct lineFile *inLf = pslFileOpen(srcPslFile);
    struct psl *psl;
    while ((psl = pslNext(inLf)) != NULL) {
        struct genebankCds *cds = getCds(psl->qName);
        struct genePred *gp = genePredFromPsl3(psl,  cds, 
                                  unsigned optFields, unsigned options,
                                  int cdsMergeSize, int utrMergeSize);
        
        pslFree(&psl);
    }

    lineFileClose(&inLf);
    return srcGenes;
}

/* create transMapInfo file */
static void augustusGenePredMk(char *cdsFile,
                               char *srcPslFile,
                               char *destPslFile,
                               char *augustusGpFile) {
    struct hash *cdsTbl = (cdsFile != NULL) ? cdsTblLoad(cdsFile) : NULL;
    struct hash *srcGenes = srcGeneLoad(cdsTbl, srcPslFile);
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 3) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    char *cdsFile = optionVal("cdsFile", NULL);
    augustusGenePredMk(cdsFile, argv[1], argv[2], argv[3]);
    return 0;
}

