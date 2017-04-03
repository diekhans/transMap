"""test configuration that can be overridden from the environment"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def getConfig(configPyFile, dataDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None, buildTmpDir=None):
    version = "V4rebuild"

    # Modify versionDir to build in different directory while using
    # version in file names. 
    versionDir = version
    
    # increment to prevent reusing batch directories on restart
    batchGen = 1
    transMapDir = os.path.join("/hive/data/outside/transMap", versionDir)
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
