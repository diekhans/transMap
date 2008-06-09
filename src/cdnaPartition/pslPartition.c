static const char *usageMsg =
    "pslPartition inPsl outPslPat\n"
    "\n"
    "split file by first n characters of accession\n"
    "outPslPat is a sprintf style format string with a single %%s\n"
    "\n"
    "Options:\n"
    "  -prefixLen=2 length off accession prefix to use\n";

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "psl.h"
#include "lib/accPartition.h"


/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"prefixLen", OPTION_INT},
    {NULL, 0},
};

/* partition psl by first two characters of accession */
static void pslPartition(int prefixLen, char *inPslFile, char *outPslPat) {
    struct accPartition *ap = accPartitionNew(prefixLen, outPslPat);
    struct lineFile *inLf = pslFileOpen(inPslFile);
    struct psl *psl;
    while ((psl = pslNext(inLf)) != NULL) {
        pslTabOut(psl, accPartitionGet(ap, psl->qName));
        pslFree(&psl);
    }
    
    accPartitionFree(ap);
    lineFileClose(&inLf);
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 3) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    pslPartition(optionInt("prefixLen", 2),
                 argv[1], argv[2]);
    return 0;
}

