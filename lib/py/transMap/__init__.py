"""standard setup for TransMap programs"""

import re


def getTwoBit(db):
    return "/hive/data/genomes/{db}/{db}.2bit".format(db=db)


def getChromSizes(db):
    return "/hive/data/genomes/{db}/chrom.sizes".format(db=db)


def parseTransMapId(transMapId):
    """parse a transmap id in the form srcDb:acc.ver-uniqmods... into a list of
    (srcDb, acc.ver, uniqmods)
    """
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+\\.[0-9]+)-([0-9.]+)$", transMapId)
    if m is None:
        raise Exception("can not parse transmap id \"{}\"".format(transMapId))
    return m.groups()


def parseAccVer(accVer):
    """parse acc.ver into (acc, ver)"""
    return accVer.split(".")
