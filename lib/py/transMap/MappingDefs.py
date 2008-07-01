from transMap import GenomeDefs
from pycbio.sys import setOps


class MappingDefs(object):
    """Get mapping definitions based on command line options"""

    @staticmethod
    def addOptions(parser):
        "add cmd line options used in creating GenomeDefs object"
        parser.add_option("--includeClade", action="append", dest="includeClades", default=None,
                          help="""Restrict to these clades,  Can be specified multiple times.""")
        parser.add_option("--otherDestDb", action="append", dest="otherDestDbs", default=None,
                          help="""Other destDbs to include beside the most current. Can be specified multiple times.""")
        parser.add_option("--excludeDb", action="append", dest="excludeDbs", default=None,
                          help="""Exclude these databases. Can be specified multiple times.""")
        parser.add_option("--onlyDestDb", action="append", dest="onlyDestDbs", default=None,
                          help="""Only generate data for these destDbs. Can be specified multiple times.""")
        parser.add_option("--cdnaType", action="append", dest="cdnaTypes", default=None,
                          help="""Only generate data for these cDNA types. Can be specified multiple times, valid values are:""" + setOps.setJoin(set(GenomeDefs.CDnaType.values), ","))

    def __init__(self, opts, inclAllMappings=False):
        "construct GenomeDefs and mappings"
        self.defs = GenomeDefs.load(opts.genomeDefsPickle)

        # restrictions; default to full set of items to make code simpler
        self.onlyDestDbs = self.defs.mapDbNamesSet(opts.onlyDestDbs) if (opts.onlyDestDbs != None) else None
        self.otherDestDbs = self.defs.mapDbNamesSet(opts.otherDestDbs)
        self.excludeDbs = self.defs.mapDbNamesSet(opts.excludeDbs)
        self.clades = setOps.mkfzset(opts.includeClades) if (opts.includeClades != None) else self.defs.clades
        if opts.cdnaTypes != None:
            # convert to Enumeration
            self.cdnaTypes = set()
            for c in opts.cdnaTypes:
                self.cdnaTypes.add(GenomeDefs.CDnaType(c))
            self.cdnaTypes = frozenset(self.cdnaTypes)
        else:
            self.cdnaTypes = self.defs.cdnaTypes

        self.mappings = self.__getMappingSets()
        self.allMappings = None
        if inclAllMappings:
            self.allMappings = self.__getMappingSets(inclMissingChains=True)

    def __useMapping(self, chainsSet, inclMissingChains):
        return ((chainsSet.srcDb not in self.excludeDbs)
                and (inclMissingChains or chainsSet.haveChains())
                and chainsSet.srcDb.matches(self.clades, self.cdnaTypes))

    def __getDestDbs(self):
        "get set of destDbs, given restrictions"
        destDbs = set(self.otherDestDbs)
        if self.onlyDestDbs != None:
            destDbs |= self.onlyDestDbs
        else:
            for org in self.defs.orgs.itervalues():
                if org.dbs[0].clade in self.clades:
                    destDbs.add(org.dbs[0])
        destDbs = destDbs.difference(self.excludeDbs)
        return destDbs

    def __getMappingSets(self, inclMissingChains=False):
        """Get a GenomeMappings object with chains describing possible
        mappings given the restrictions.  Only the latest db for an organism
        is used unless an earlier db is specified in otherDestDbs"""
        gms = GenomeMappings()
        for destDb in self.__getDestDbs():
            for ent in destDb.srcChainsSets.iterentries():
                # ent[0] is newest srcDb chainsSet
                if self.__useMapping(ent[0], inclMissingChains):
                    gms._addMapping(ent[0])
        gms.finish()
        return gms


class GenomeMappings(list):
    "list of ChainsSet objects, along with sets of databases, to map"
    def __init__(self):
        self.srcDbs = set()
        self.destDbs = set()
        self.chainsMap = set()  # unordered, but quick lookup, self is sorted

    def _addMapping(self, chains):
        "add a mapping"
        self.append(chains)
        self.chainsMap.add(chains)
        self.srcDbs.add(chains.srcDb)
        self.destDbs.add(chains.destDb)

    @staticmethod
    def __srcDestCmp(a, b):
        """sort comparison by srcDb and destDb fields objects"""
        d = cmp(a.srcDb.orgDbName, b.srcDb.orgDbName)
        if d == 0:
            d = cmp(a.srcDb.dbNum, b.srcDb.dbNum)
            if d == 0:
                d = cmp(a.destDb.orgDbName, b.destDb.orgDbName)
                if d == 0:
                    d = cmp(a.destDb.dbNum, b.destDb.dbNum)
        return d


    def finish(self):
        "sort by srcDb, destDb"
        self.sort(self.__srcDestCmp)

    def dump(self, fh):
        "print for debugging purposes"
        prRowv(fh, "Mappings:")
        for ch in self:
            fh.write("\t")
            fh.write(str(ch))
            fh.write("\n")

