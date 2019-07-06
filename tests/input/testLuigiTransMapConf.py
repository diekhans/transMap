"""test configuration that can be overridden from the environment"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def getConfig(configPyFile, dataRootDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None):
    if dataRootDir is None:
        dataRootDir = "luigiTest/root"
    conf = TransMapConf(configPyFile,
                        dataRootDir=dataRootDir,
                        srcHgDb=srcHgDb,
                        destHgDb=destHgDb,
                        annotationType=annotationType,
                        chainType=chainType,
                        version="V7",
                        batchGen=1)
    return conf
