#ifndef alignPair_h
#define alignPair_h
struct lineFile;

/* pairing of psl and mapInfo */
struct alignPair {
    struct alignPair *next;
    struct psl *psl;
    struct mapInfo *mi;
    char *key;    // lazy construction my alignPairGetKey
};

/* constructor */
struct alignPair *alignPairNew(struct psl *psl, struct mapInfo *mi);

/* destructor */
void alignPairFree(struct alignPair **apPtr);

/* construct a key for an alignment pair */
char* alignPairMkKey(char *qName, int qStart, int qEnd,
                     char *tName, int tStart, int tEnd);

/* construct a key for an alignment pair from a psl */
char* alignPairMkPslKey(struct psl *psl);

/* construct a key for an alignment pair from a mapInfo object */
char* alignPairMkMapInfoKey(struct mapInfo *mi);

/* get key, creating on first use */
INLINE char *alignPairGetKey(struct alignPair *ap) {
    if (ap->key == NULL) {
        ap->key = alignPairMkPslKey(ap->psl);
    }
    return ap->key;
}

/* compare key fields in two alignPair objects, also sorting my query order */
boolean alignPairKeysCmp(struct alignPair *a, struct alignPair *b);

/* load an alignPair object from two parallel lineFiles */
struct alignPair *alignPairLoad(struct lineFile *pslLf, struct lineFile *miLf);

/* load all alignPair objects from two parallel files */
struct alignPair *alignPairLoadAll(char *pslFile, char *miFile);

/* write an alignPair object to two parallel files */
void alignPairWrite(struct alignPair *ap, FILE *pslFh, FILE *miFh);

/* write all alignPair objects to two parallel files */
void alignPairWriteAll(struct alignPair *aps, char *pslFile, char *miFile);

/* sort align pair objects, first by query, then by target; this is also key order */
int alignPairQuerySortCmp(const void *va, const void *vb);

#endif 
