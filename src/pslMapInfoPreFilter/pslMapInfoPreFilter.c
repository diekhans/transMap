static const char *usageMsg =
    "pslMapInfoPreFilter [options] pslIn mapInfoIn pslOut mapInfoOut\n"
    "\n"
    "Prefiltering of psl and mapInfo files immediately after pslMap.  This deals with\n"
    "ensuring that psls and mapInfo are unique based on the tuple:\n"
    "    qName, qStart, qEnd, tName, tStart, tEnd, strand\n"
    "This allows down-stream filter to use pslMapInfoSelect and get rid of duplicates based\n"
    "on chains that a identical in the gene region\n"
    "\n"
    "Options:\n"
    "   -chainSubset=val - set chain subset to this value\n";

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "localmem.h"
#include "psl.h"
#include "options.h"
#include "lib/mapInfo.h"
#include "lib/alignPair.h"

/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {"chainSubset", OPTION_STRING},
    {NULL, 0}
};
static boolean haveChainSubset = FALSE;
static enum mapInfoChainSubset chainSubset = mapInfoUnknown;

/* get an entry from the key-sorted list; if there are multiple with the same key,
 * keep return the best scoring and drop the rest */
static struct alignPair *popBestNextKey(struct alignPair **apsIn) {
    // don't worry about leaking memory
    if (*apsIn == NULL) {
        return NULL;
    }
    struct alignPair *bestAp = slPopHead(apsIn);
    int bestScore = pslScore(bestAp->psl);
    while ((*apsIn != NULL) &&  (alignPairKeysCmp(bestAp, *apsIn) == 0)) {
        struct alignPair *ap = slPopHead(apsIn);
        int score = pslScore(ap->psl);
        if (score > bestScore) {
            bestAp = ap;
            bestScore = score;
        }
    }
    return bestAp;
}

/* pre-filter psl/mapInfo files */
static void pslMapInfoPreFilter(char *pslIn,
                                char *mapInfoIn,
                                char *pslOut,
                                char *mapInfoOut) {
    struct alignPair *apsIn = alignPairLoadAll(pslIn, mapInfoIn);
    slSort(&apsIn, alignPairQuerySortCmp);  // key order!!
    
    struct alignPair *apsOut = NULL, *ap;
    while ((ap = popBestNextKey(&apsIn)) != NULL) {
        slAddHead(&apsOut, ap);
        if (haveChainSubset) {
            ap->mi->chainSubset = chainSubset;
        }
    }
    slReverse(&apsOut);
    alignPairWriteAll(apsOut, pslOut, mapInfoOut);    
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    haveChainSubset = optionExists("chainSubset");
    if (haveChainSubset) {
        chainSubset =  mapInfoParseChainSubset(optionVal("chainSubset", NULL));
    }
    pslMapInfoPreFilter(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}

