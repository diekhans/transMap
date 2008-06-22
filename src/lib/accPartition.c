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

/* get path for a partition. WARNING: static return */
char *accPartitionPath(struct accPartition *ap,
                       char* accPre) {
    static char path[PATH_LEN];
    safef(path, sizeof(path), ap->pathPat, accPre);
    return path;
}

/* obtain file associate with partition, opening if needed */
FILE *accPartitionGet(struct accPartition *ap,
                      char* acc) {
    char accPre[ap->prefixLen+1];
    strncpy(accPre, acc, ap->prefixLen);
    accPre[ap->prefixLen] = '\0';
    struct hashEl *hel = hashStore(ap->files, accPre);
    if (hel->val == NULL) {
        hel->val = mustOpen(accPartitionPath(ap, accPre), "w");
    }
    return hel->val;
}

/* obtain list of prefixes */
struct slName *accPartitionPrefixes(struct accPartition *ap) {
    struct slName *pres = NULL;
    struct hashCookie cookie = hashFirst(ap->files);
    struct hashEl *hel;
    while ((hel = hashNext(&cookie)) != NULL) {
        slSafeAddHead(&pres, slNameNew(hel->name));
    }
    slNameSort(&pres);
    return pres;
}
