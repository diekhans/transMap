from transMap import GenomeDefs


class MappingDefs(object):
    """Get mapping definitions based on command line options"""

    @staticmethod
    def addOptions(parser):
        "add cmd line options used in creating GenomeDefs object"
        parser.add_option("--clade", action="append", dest="clades", default=None,
                          help="""Restrict to these clades,  Can be specified multiple times.""")
        parser.add_option("--otherDestDb", action="append", dest="otherDestDbs", default=None,
                          help="""Other destDbs to include beside the most current. Can be specified multiple times.""")
        parser.add_option("--skipDb", action="append", dest="skipDbs", default=None,
                          help="""Skip these databases. Can be specified multiple times.""")
        parser.add_option("--onlyDestDb", action="append", dest="onlyDestDbs", default=None,
                          help="""Only generate data for these destDbs. Can be specified multiple times.""")

    def __init__(self, opts, inclAllMappings=False):
        "construct GenomeDefs and mappings"
        self.defs = GenomeDefs.load(opts.genomeDefsPickle)
        self.mappings = self.defs.getMappingSets(onlyDestDbs=opts.onlyDestDbs,
                                            otherDestDbs=opts.otherDestDbs, 
                                            clades=opts.clades,
                                            skipDbs=opts.skipDbs)
        self.allMappings = None
        if inclAllMappings:
            self.allMappings = self.defs.getMappingSets(onlyDestDbs=opts.onlyDestDbs,
                                                        otherDestDbs=opts.otherDestDbs, 
                                                        clades=opts.clades,
                                                        skipDbs=opts.skipDbs,
                                                        inclMissingChains=True,
                                                        inclRevChains=True)
