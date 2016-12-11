"""test configuration that can be overridden from the environment"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def getConfig(configPyFile, dataDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None, buildTmpDir=None):
    if dataDir is None:
        dataDir = "luigiTest/data"
    if buildTmpDir is None:
        buildTmpDir = "luigiTest/build"
    conf = TransMapConf(configPyFile,
                        dataDir=dataDir,
                        srcHgDb=srcHgDb,
                        destHgDb=destHgDb,
                        annotationType=annotationType,
                        chainType=chainType,
                        buildTmpDir=buildTmpDir,
                        version="V7")
    conf.requiredPreviousDestHgDbs = frozenset()
    return conf
