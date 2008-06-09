#include "common.h"
#include "lib/accPartition.h"
#include "hash.h"

/* construct a new object */
struct accPartition *accPartitionNew(int prefixLen,
                                     char *pathPat) {
    struct accPartition *ap;
    AllocVar(ap);
    ap->prefixLen = prefixLen;
    safecpy(ap->pathPat, sizeof(ap->pathPat), pathPat);
    ap->files = hashNew(15);
    return ap;
}

/* free object, closing files */
void accPartitionFree(struct accPartition *ap) {
    if (ap != NULL) {
        struct hashCookie cookie = hashFirst(ap->files);
        struct hashEl *hel;
        while ((hel = hashNext(&cookie)) != NULL) {
            FILE *fh = hel->val;
            carefulClose(&fh);
        }
        hashFree(&ap->files);
        freeMem(ap);
    }
}

/* obtain file associate with partition, opening if needed */
FILE *accPartitionGet(struct accPartition *ap,
                      char* acc) {
    char accPre[ap->prefixLen+1];
    strncpy(accPre, acc, ap->prefixLen);
    accPre[ap->prefixLen] = '\0';
    struct hashEl *hel = hashStore(ap->files, accPre);
    if (hel->val == NULL) {
        char path[PATH_LEN];
        safef(path, sizeof(path), ap->pathPat, accPre);
        hel->val = mustOpen(path, "w");
    }
    return hel->val;
}
