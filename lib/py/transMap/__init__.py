"""standard setup for TransMap programs"""

import re
import os
import sys


def setSortLocale():
    "make sure sort command does the right thing"
    os.environ["LC_ALL"] = "C"


def getSortProg():
    "requires gnu sort"
    return "gsort" if sys.platform == 'darwin' else "sort"


def getTwoBit(db):
    paths = ["/scratch/data/{db}/{db}.2bit".format(db=db),
             "/hive/data/genomes/{db}/{db}.2bit".format(db=db)]
    for twoBit in paths:
        if os.path.exists(twoBit):
            return twoBit
    raise Exception("can't find 2bit for {db}".format(db=db))


def getChromSizes(db):
    return "/hive/data/genomes/{db}/chrom.sizes".format(db=db)


def parseAlignId(alignId):
    """parse a transmap id in the form srcDb:acc.ver-uniqmods... into a list of
    (srcDb, acc.ver, uniqmods).  This will parse a srouce or target alignment id
    """
    # older versions of ensembl tracks didn't have versions
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+(\\.[0-9]+)?)-([0-9.]+)$", alignId)
    if m is None:
        raise ValueError("can not parse transmap alignment id \"{}\"".format(alignId))
    return (m.group(1), m.group(2), m.group(4))


def alignIdDropOneUniq(alignId):
    "drop one level of alignment id uniqueness"
    parts = parseAlignId(alignId)
    idot = parts[2].rfind('.')
    if idot < 0:
        return "{}:{}".format(parts[0], parts[1])
    else:
        return "{}:{}-{}".format(parts[0], parts[1], parts[2][0:idot])


def alignIdToSrcId(alignId):
    "remove the unique suffix only, return srcId, including srcDb"
    # older versions of ensembl tracks didn't have versions
    m = re.match("^([a-zA-Z0-9]+:[A-Z0-9_]+(\\.[0-9]+)?)-([0-9.]+)$", alignId)
    if m is None:
        raise ValueError("can not parse alignment id \"{}\"".format(alignId))
    return m.group(1)


def srcIdToAccv(srcId):
    "convert a srcId to an accv"
    # older versions of ensembl tracks didn't have versions
    m = re.match("^([a-zA-Z0-9]+):([A-Z0-9_]+(\\.[0-9]+)?)$", srcId)
    if m is None:
        raise ValueError("can not parse source id \"{}\"".format(srcId))
    return m.group(2)


def parseAccVer(accver):
    """parse acc.ver into (acc, ver)"""
    m = re.match("^([A-Z0-9_]+)\\.([0-9]+)$", accver)
    if m is None:
        raise ValueError("can not parse accession.version id \"{}\"".format(accver))
    return (m.group(1), m.group(2))


def accvToAcc(accver):
    """parse acc.ver into (acc, ver)"""
    return parseAccVer(accver)[0]


_hgDbNameRe = re.compile("^([a-zA-Z]+)([0-9]+)$")


def hgDbNameParse(hgDb):
    """parse a hg database name into (orgDbName, orgDbNum)"""
    m = _hgDbNameRe.match(hgDb)
    if m is None:
        raise Exception("can't parse UCSC database name: {}".format(hgDb))
    return (m.group(1), int(m.group(2)))


def hgDbNameCheck(hgDb):
    """check a hg database name is a valid name"""
    return _hgDbNameRe.match(hgDb) is not None
