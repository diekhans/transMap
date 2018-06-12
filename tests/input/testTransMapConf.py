"""test configuration that can be overridden from the environment"""

import os
from transMap.transMapConf import TransMapConf
from transMap.genomeData import ChainType


def _getConfVal(name, supplied):
    """return configuration value if supplied is not None
    else try environment"""
    if supplied is not None:
        return supplied
    else:
        return os.environ.get(name, None)

def getConfig(configPyFile, dataRootDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None):
    conf = TransMapConf(configPyFile,
                        dataRootDir=_getConfVal("dataRootDir", dataRootDir),
                        srcHgDb=_getConfVal("srcHgDb", srcHgDb),
                        destHgDb=_getConfVal("destHgDb", destHgDb),
                        annotationType=_getConfVal("annotationType", annotationType),
                        chainType=_getConfVal("chainType", chainType),
                        version="V1",
                        batchGen=1)
    if os.environ.get("noRequired", None) is not None:
        conf.requiredPreviousDestHgDbs = frozenset()
    conf.gbdbDir = _getConfVal("gbdbDir", None)   # None will prevent overwrite of real /gbdb
    return conf
