static const char *usageMsg =
    "pslPartition inPsl outPslPat\n"
    "\n"
    "split file by first n characters of accession\n"
    "outPslPat is a sprintf style format string with a single %%s\n"
    "\n"
    "Options:\n"
    "  -prefixLen=2 length off accession prefix to use\n"
    "  -prefixList=out - write list of prefix to this file.\n"
    "  -pslList=out - write list of psls create to this file.\n";

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "psl.h"
#include "hash.h"
#include "lib/accPartition.h"


/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"prefixLen", OPTION_INT},
    {"prefixList", OPTION_STRING},
    {"pslList", OPTION_STRING},
    {NULL, 0},
};

/* partition PSLs */
static void partitionPslFile(int prefixLen, char *inPslFile,
                             struct accPartition *ap) {
    struct lineFile *inLf = pslFileOpen(inPslFile);
    struct psl *psl;
    while ((psl = pslNext(inLf)) != NULL) {
        pslTabOut(psl, accPartitionGet(ap, psl->qName));
        pslFree(&psl);
    }
    lineFileClose(&inLf);
}

/* write PSL list */
static void listPsls(struct accPartition *ap,
                    char *pslList) {
    FILE *fh = mustOpen(pslList, "w");
    struct slName *pres = accPartitionPrefixes(ap);
    for (struct slName *pre = pres; pre != NULL; pre = pre->next) {
        fprintf(fh, "%s\n", accPartitionPath(ap, pre->name));
    }
    slFreeList(&pres);
    carefulClose(&fh);
}

/* write prefixes list */
static void listPrefixes(struct accPartition *ap,
                         char *prefixList) {
    FILE *fh = mustOpen(prefixList, "w");
    struct slName *pres = accPartitionPrefixes(ap);
    for (struct slName *pre = pres; pre != NULL; pre = pre->next) {
        fprintf(fh, "%s\n", pre->name);
    }
    slFreeList(&pres);
    carefulClose(&fh);
}

/* partition psl by first two characters of accession */
static void pslPartition(int prefixLen, char *inPslFile, char *outPslPat,
                         char *prefixList, char *pslList) {
    struct accPartition *ap = accPartitionNew(prefixLen, outPslPat);
    partitionPslFile(prefixLen, inPslFile, ap);
    if (prefixList != NULL) {
        listPrefixes(ap, prefixList);
    }
    if (pslList != NULL) {
        listPsls(ap, pslList);
    }
    accPartitionFree(ap);
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 3) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    pslPartition(optionInt("prefixLen", 2), argv[1], argv[2],
                 optionVal("prefixList", NULL),
                 optionVal("pslList", NULL));
    return 0;
}

