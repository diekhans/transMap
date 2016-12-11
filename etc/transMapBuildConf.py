"""test configuration that can be overridden from the environment"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def getConfig(configPyFile, dataDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None, buildTmpDir=None):
    version = "V4"
    # increment to prevent reusing batch directories on restart
    batchGen = 9
    transMapDir = os.path.join("/hive/data/outside/transMap", version)
    dataDir = os.path.join(transMapDir, "data")
    buildTmpDir = os.path.join(transMapDir, "build.tmp")
    conf = TransMapConf(configPyFile,
                        dataDir=dataDir,
                        srcHgDb=srcHgDb,
                        destHgDb=destHgDb,
                        annotationType=annotationType,
                        chainType=chainType,
                        buildTmpDir=buildTmpDir,
                        version=version,
                        batchGen=batchGen)
    return conf
