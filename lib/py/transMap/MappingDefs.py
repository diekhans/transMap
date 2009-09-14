from transMap import GenomeDefs

from pycbio.sys import setOps,fileOps
from pycbio.sys.MultiDict import MultiDict

class MappingDefsBuild(object):
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

    def __init__(self, opts):
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
        else:
            self.cdnaTypes = self.defs.cdnaTypes
        self.cdnaTypes = frozenset(self.cdnaTypes)

        self.mappings = self.__getMappingSets()

    def __inclDestOrg(self, org):
        "should a destDb be used from this org?"
        return org.dbs[0].clade in self.clades

    def __getNewestDestDb(self, org):
        "get the newest, skipping excluded"
        # ordered newest to oldest
        for db in org.dbs:
            if db not in self.excludeDbs:
                return db
        return None

    def __getDestDbs(self):
        "get set of destDbs, given restrictions"
        destDbs = set(self.otherDestDbs)
        if self.onlyDestDbs != None:
            destDbs |= self.onlyDestDbs
        else:
            for org in self.defs.orgs.itervalues():
                if self.__inclDestOrg(org):
                    db = self.__getNewestDestDb(org)
                    if db != None:
                        destDbs.add(db)
        return destDbs

    def __useMapping(self, chainsSet):
        return ((chainsSet.srcDb not in self.excludeDbs) and chainsSet.haveChains()
                and chainsSet.srcDb.matches(self.clades, self.cdnaTypes))

    def __getSrcOrgDestDbMappings(self, srcOrg, destDb, cdnaType, chainsSets, gms):
        # chainsSets are sorted new to oldest srcDb
        for chainsSet in chainsSets:
            if self.__useMapping(chainsSet) and (cdnaType in chainsSet.srcDb.cdnaTypes):
                gms.addMapping(chainsSet, cdnaType)
                break # stop a newests

    def __getDestDbMappings(self, destDb, gms):
        "build all GenomeMapping objects for a destDb"
        for srcOrg in destDb.srcChainsSets.iterkeys():
            for cdnaType in GenomeDefs.CDnaType.values:
                self.__getSrcOrgDestDbMappings(srcOrg, destDb, cdnaType, destDb.srcChainsSets[srcOrg], gms)

    def __getMappingSets(self):
        """Get a GenomeMappings object with chains describing possible
        mappings given the restrictions."""
        gms = GenomeMappings()
        for destDb in self.__getDestDbs():
            self.__getDestDbMappings(destDb, gms)
        return gms

class GenomeMapping(object):
    "mapping of a cDNA type from a source to desc db"
    def __init__(self, chainsSet, cdnaType):
        self.srcDb = chainsSet.srcDb
        self.cdnaType = cdnaType
        self.destDb = chainsSet.destDb
        self.chainsSet = chainsSet

class DestDbMappings(object):
    "mappings for one destDb, split by cDNA type"
    def __init__(self, destDb):
        self.destDb = destDb
        self.srcDbs = set()
        self.mappingsBySrcDb = MultiDict()
        
    def addMapping(self, gm):
        "add a genome mapping"
        assert(gm.destDb == self.destDb)
        self.srcDbs.add(gm.srcDb)
        self.mappingsBySrcDb.add(gm.srcDb, gm)

    def __printSrcDbMappings(self, fh, srcDb):
        gms = self.mappingsBySrcDb[srcDb]
        cdnaTypes = frozenset([gm.cdnaType for gm in gms])
        chainTypes = frozenset(gms[0].chainsSet.byType.iterkeys())
        fileOps.prRowv(fh, self.destDb.name, srcDb.name, setOps.setJoin(cdnaTypes, ","),
                       setOps.setJoin(chainTypes, ","), gms[0].chainsSet.dist)

    def printMappings(self, fh):
        srcDbs = list(self.srcDbs)
        srcDbs.sort(key=lambda a: a.name)
        for srcDb in srcDbs:
            self.__printSrcDbMappings(fh, srcDb)
        

class GenomeMappings(list):
    "list of DestDbMappings objects for all selected destination genomes"
    def __init__(self):
        self.destDbMappingsTbl = dict()
        self.bySrcDb = MultiDict()

    def __obtainDestDb(self, destDb):
        destDbMappings = self.destDbMappingsTbl.get(destDb)
        if destDbMappings == None:
            destDbMappings = self.destDbMappingsTbl[destDb] =  DestDbMappings(destDb)
            self.append(destDbMappings)
        return destDbMappings
    
    def addMapping(self, chainsSet, cdnaType):
        "add a mapping"
        destDbMappings = self.__obtainDestDb(chainsSet.destDb)
        gm = GenomeMapping(chainsSet, cdnaType)
        destDbMappings.addMapping(gm)
        self.bySrcDb.add(gm.srcDb, gm)

    def printMappings(self, fh):
        destDbs = list(self.destDbMappingsTbl.iterkeys())
        destDbs.sort(key=lambda a: a.name)
        for destDb in destDbs:
            self.destDbMappingsTbl[destDb].printMappings(fh)
