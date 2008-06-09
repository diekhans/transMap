static const char *usageMsg =
    "faPartition inFa outFaPat\n"
    "\n"
    "split file by first n characters of accession\n"
    "outFaPat is a sprintf style format string with a single %%s\n"
    "\n"
    "Options:\n"
    "  -prefixLen=2 length off accession prefix to use\n";

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "fa.h"
#include "lib/accPartition.h"


/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"prefixLen", OPTION_INT},
    {NULL, 0},
};

/* partition one sequence */
static void faPartitionSeq(struct accPartition *ap, DNA *seq, int size, char *header) {
    char buf[256];
    safecpy(buf, sizeof(buf), header);
    char *acc = firstWordInLine(buf);
    faWriteNext(accPartitionGet(ap, acc), header, seq, size);
}

/* partition fasta by first two characters of accession */
static void faPartition(int prefixLen, char *inFaFile, char *outFaPat) {
    struct accPartition *ap = accPartitionNew(prefixLen, outFaPat);
    struct lineFile *inLf = lineFileOpen(inFaFile, TRUE);
    DNA *seq;
    int size;
    char *header;
    while (faSpeedReadNext(inLf, &seq, &size, &header)) {
        faPartitionSeq(ap, seq, size, header);
    }
    
    accPartitionFree(ap);
    faFreeFastBuf();
    lineFileClose(&inLf);
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 3) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    faPartition(optionInt("prefixLen", 2),
                 argv[1], argv[2]);
    return 0;
}

