/* mapInfo.h was originally generated by the autoSql program, which also 
 * generated mapInfo.c and mapInfo.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef MAPINFO_H
#define MAPINFO_H

#define MAPINFO_NUM_COLS 27

struct mapInfo
/* psl mapping information */
    {
    struct mapInfo *next;  /* Next in singly linked list. */
    char *srcQName;	/* srcQName */
    int srcQStart;	/* srcQStart */
    int srcQEnd;	/* srcQEnd */
    int srcQSize;	/* srcQSize */
    char *srcTName;	/* srcTName */
    int srcTStart;	/* srcTStart */
    int srcTEnd;	/* srcTEnd */
    char srcStrand;	/* srcStrand */
    int srcAligned;	/* srcAligned */
    char *mappingQName;	/* mappingQName */
    int mappingQStart;	/* mappingQStart */
    int mappingQEnd;	/* mappingQEnd */
    char *mappingTName;	/* mappingTName */
    int mappingTStart;	/* mappingTStart */
    int mappingTEnd;	/* mappingTEnd */
    char mappingStrand;	/* mappingStrand */
    char *mappingId;	/* mappingId */
    char *mappedQName;	/* mappedQName */
    int mappedQStart;	/* mappedQStart */
    int mappedQEnd;	/* mappedQEnd */
    char *mappedTName;	/* mappedTName */
    int mappedTStart;	/* mappedTStart */
    int mappedTEnd;	/* mappedTEnd */
    char mappedStrand;	/* mappedStrand */
    int mappedAligned;	/* mappedAligned */
    int qStartTrunc;	/* qStartTrunc */
    int qEndTrunc;	/* qEndTrunc */
    };

void mapInfoStaticLoad(char **row, struct mapInfo *ret);
/* Load a row from mapInfo table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct mapInfo *mapInfoLoad(char **row);
/* Load a mapInfo from row fetched with select * from mapInfo
 * from database.  Dispose of this with mapInfoFree(). */

struct mapInfo *mapInfoLoadAll(char *fileName);
/* Load all mapInfo from whitespace-separated file.
 * Dispose of this with mapInfoFreeList(). */

struct mapInfo *mapInfoLoadAllByChar(char *fileName, char chopper);
/* Load all mapInfo from chopper separated file.
 * Dispose of this with mapInfoFreeList(). */

#define mapInfoLoadAllByTab(a) mapInfoLoadAllByChar(a, '\t');
/* Load all mapInfo from tab separated file.
 * Dispose of this with mapInfoFreeList(). */

struct mapInfo *mapInfoCommaIn(char **pS, struct mapInfo *ret);
/* Create a mapInfo out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new mapInfo */

void mapInfoFree(struct mapInfo **pEl);
/* Free a single dynamically allocated mapInfo such as created
 * with mapInfoLoad(). */

void mapInfoFreeList(struct mapInfo **pList);
/* Free a list of dynamically allocated mapInfo's */

void mapInfoOutput(struct mapInfo *el, FILE *f, char sep, char lastSep);
/* Print out mapInfo.  Separate fields with sep. Follow last field with lastSep. */

#define mapInfoTabOut(el,f) mapInfoOutput(el,f,'\t','\n');
/* Print out mapInfo as a line in a tab-separated file. */

#define mapInfoCommaOut(el,f) mapInfoOutput(el,f,',',',');
/* Print out mapInfo as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* MAPINFO_H */
