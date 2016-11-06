"""standard setup for TransMap programs"""

import re


def getTwoBit(db):
    return "/hive/data/genomes/{db}/{db}.2bit".format(db=db)


def getChromSizes(db):
    return "/hive/data/genomes/{db}/chrom.sizes".format(db=db)


def parseAlignId(alignId):
    """parse a transmap id in the form srcDb:acc.ver-uniqmods... into a list of
    (srcDb, acc.ver, uniqmods).  This will parse a srouce or target alignment id
    """
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+\\.[0-9]+)-([0-9.]+)$", alignId)
    if m is None:
        raise Exception("can not parse transmap id \"{}\"".format(alignId))
    return m.groups()


def alignIdToSrcId(alignId):
    "remove the unique suffix only, return srcId, including srcDb"
    m = re.match("^([a-zA-Z0-9]+:[A-Z0-9_]+\\.[0-9]+)-([0-9.]+)$", alignId)
    if m is None:
        raise Exception("can not parse alignment id \"{}\"".format(alignId))
    return m.group(1)


def srcIdToAccv(srcId):
    "convert a srcId to an accv"
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+\\.[0-9]+)$", srcId)
    if m is None:
        raise Exception("can not parse source id \"{}\"".format(srcId))
    return m.group(2)


def parseAccVer(accver):
    """parse acc.ver into (acc, ver)"""
    return accver.split(".")
