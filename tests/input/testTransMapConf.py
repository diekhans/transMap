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

def getConfig(configPyFile, dataDir=None, srcHgDb=None, destHgDb=None,
              annotationType=None, chainType=None, buildTmpDir=None):
    conf = TransMapConf(configPyFile,
                        dataDir=_getConfVal("dataDir", dataDir),
                        srcHgDb=_getConfVal("srcHgDb", srcHgDb),
                        destHgDb=_getConfVal("destHgDb", destHgDb),
                        annotationType=_getConfVal("annotationType", annotationType),
                        chainType=_getConfVal("chainType", chainType),
                        buildTmpDir=_getConfVal("buildTmpDir", buildTmpDir))


    # test environ doesn't have syn chains
    conf.cladePreferedChains[frozenset(["mammal", "mammal"])] =  (ChainType.syn,  ChainType.all)
    return conf
