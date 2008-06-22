#ifndef accPartition_h
#define accPartition_h

/* object partition by acc prefix, keeps hash of open files */
struct accPartition {
    int prefixLen;
    char pathPat[PATH_LEN]; // contains a single %s to substitute with prefix
    struct hash *files;    // prefix to open FILE
};

/* construct a new object */
struct accPartition *accPartitionNew(int prefixLen,
                                     char *pathPat);

/* free object, closing files */
void accPartitionFree(struct accPartition *ap);

/* get path for a partition. WARNING: static return */
char *accPartitionPath(struct accPartition *ap,
                       char* accPre);

/* obtain file associate with partition, opening if needed */
FILE *accPartitionGet(struct accPartition *ap,
                      char* acc);

/* obtain list of prefixes */
struct slName *accPartitionPrefixes(struct accPartition *ap);

#endif
