

class TransMapConf(object):
    """Configuration class for transmap. config file should define
    a function:
       getConfig(srcHgDb=None, destHgDb=None)
    The srcHgDb  and destHgDb parameters are specified only in commands
    operating
    """
    def __init__(self, srcHgDb=None, destHgDb=None):
        self.srcHgDb = srcHgDb
        self.destHgDb = destHgDb

        # arguments to pslCDnaFilter
        self.mappedMinQSize = 20
        self.mappedMinCover = 0.20
        self.mappedMaxAligns = 20
        self.mappedGlobalNearBest = 0.01
