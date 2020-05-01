"""configuration file for transMap batch run"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def getConfig(configPyFile, dataRootDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None):
    version = "V5"
    prevVersion = "V4"

    # increment to prevent reusing batch directories on restart
    batchGen = 1
    dataRootDir = os.path.join("/hive/data/inside/transMap", version)
    prevDataRootDir = os.path.join("/hive/data/inside/transMap", prevVersion)
    conf = TransMapConf(configPyFile,
                        dataRootDir=dataRootDir,
                        srcHgDb=srcHgDb,
                        destHgDb=destHgDb,
                        annotationType=annotationType,
                        chainType=chainType,
                        version=version,
                        batchGen=batchGen,
                        prevVersion=prevVersion,
                        prevDataRootDir=prevDataRootDir)
    return conf
