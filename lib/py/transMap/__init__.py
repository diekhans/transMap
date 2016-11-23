"""standard setup for TransMap programs"""

import re
import os


def getTwoBit(db):
    paths = ["/scratch/data/{db}/{db}.2bit".format(db=db),
             "/hive/data/genomes/{db}/{db}.2bit".format(db=db)]
    for twoBit in paths:
        if os.path.exists(twoBit):
            return twoBit
    raise Exception("can't find 2bit for db{}".format(db=db))


def getChromSizes(db):
    return "/hive/data/genomes/{db}/chrom.sizes".format(db=db)


def parseAlignId(alignId):
    """parse a transmap id in the form srcDb:acc.ver-uniqmods... into a list of
    (srcDb, acc.ver, uniqmods).  This will parse a srouce or target alignment id
    """
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+\\.[0-9]+)-([0-9.]+)$", alignId)
    if m is None:
        raise Exception("can not parse transmap alignment id \"{}\"".format(alignId))
    return m.groups()

def alignIdDropOneUniq(alignId):
    "drop one level of alignment id uniqueness"
    parts = parseAlignId(alignId)
    idot = parts[1].rindex('.')
    if idot < 0:
        return parts[0]
    else:
        return "{}-{}".format(parts[0], parts[1][0:idot])


def alignIdToSrcId(alignId):
    "remove the unique suffix only, return srcId, including srcDb"
    return parseAlignId(alignId)[0]


def srcIdToAccv(srcId):
    "convert a srcId to an accv"
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+\\.[0-9]+)$", srcId)
    if m is None:
        raise Exception("can not parse source id \"{}\"".format(srcId))
    return m.group(2)


def parseAccVer(accver):
    """parse acc.ver into (acc, ver)"""
    return accver.split(".")


def accvToAcc(accver):
    """parse acc.ver into (acc, ver)"""
    return parseAccVer(accver)[0]


def hgDbNameParse(hgDb):
    """parse a hg database name into (orgDbName, orgDbNum)"""
    m = re.match("^([a-zA-Z]+)([0-9]+)$", hgDb)
    if m is None:
        raise Exception("can't parse UCSC database name: " + hgDb)
    return (m.group(1), int(m.group(2)))
