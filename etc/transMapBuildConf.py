"""test configuration that can be overridden from the environment"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def getConfig(configPyFile, dataRootDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None):
    version = "V5"

    # Modify versionDir to build in different directory while using
    # version in file names.
    versionDir = version

    # increment to prevent reusing batch directories on restart
    batchGen = 1
    dataRootDir = os.path.join("/hive/data/inside/transMap", versionDir)
    conf = TransMapConf(configPyFile,
                        dataRootDir=dataRootDir,
                        srcHgDb=srcHgDb,
                        destHgDb=destHgDb,
                        annotationType=annotationType,
                        chainType=chainType,
                        version=version,
                        batchGen=batchGen)
    return conf
