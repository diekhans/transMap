static const char *usageMsg =
    "transMapInfoMk transMapInfo mapInfo1 [mapInfo2 ...]\n"
    "\n"
    "Generate transMapInfo file from a set of mapInfo files. Assumes\n"
    "that directory name of each mapInfo is source genome database.\n";

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "localmem.h"
#include "psl.h"
#include "options.h"
#include "transMapInfo.h"
#include "lib/mapInfo.h"

/* command line options and values */
static struct optionSpec optionSpecs[] = {
    {NULL, 0}
};

/* parse data from a path: WARNING: static return */
static char *parseDbFromPath(char *path) {
    static char buf[PATH_LEN];
    safecpy(buf, sizeof(buf), path);
    char *slash = strrchr(buf, '/');
    if (slash != NULL) {
        *slash = '\0';
        slash = strrchr(buf, '/');
        if (slash != NULL) {
            return slash+1;
        }
    }
    errAbort("can't parse DB from path: %s", path);
    return NULL;
}

/* convert and output one mapInfo object */
static void processMapInfo(FILE *outFh,
                           char *srcDb,
                           struct mapInfo *mi) {
    struct transMapInfo tmi;
    ZeroVar(&tmi);
    safecpy(tmi.srcDb, sizeof(tmi.srcDb), srcDb);
    tmi.srcId = mi->srcQName;
    tmi.mappingId = mi->mappingId;
    tmi.mappedId = mi->mappedQName;
    transMapInfoTabOut(&tmi, outFh);
}

/* process one mapInfo file */
static void mapInfoToTransMapInfo(FILE *outFh,
                                  char *mapInfoIn) {
    char *srcDb = parseDbFromPath(mapInfoIn);
    struct lineFile *inLf = lineFileOpen(mapInfoIn, TRUE);
    char *row[MAPINFO_NUM_COLS];
    while (lineFileNextRowTab(inLf, row, MAPINFO_NUM_COLS)) {
        struct mapInfo *mi = mapInfoLoad(row);
        processMapInfo(outFh, srcDb, mi);
        mapInfoFree(&mi);
    }
    lineFileClose(&inLf);
}

/* create transMapInfo file */
static void transMapInfoMk(char *transMapInfoOut,
                           int mapInfoCnt,
                           char **mapInfoIns) {
    int i;
    FILE *outFh = mustOpen(transMapInfoOut, "w");
    for (i = 0; i < mapInfoCnt; i++) {
        mapInfoToTransMapInfo(outFh, mapInfoIns[i]);
    }
    carefulClose(&outFh);
}

/* entry */
int main(int argc, char **argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc < 3) {
        errAbort("wrong # of args: %s", usageMsg);
    }
    transMapInfoMk(argv[1], argc-2, argv+2);
    return 0;
}

