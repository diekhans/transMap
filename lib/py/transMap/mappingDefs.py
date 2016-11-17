from . import genomeData

from pycbio.sys import setOps,fileOps
from pycbio.sys.multiDict import MultiDict

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
        parser.add_option("--onlySrcDb", action="append", dest="onlySrcDbs", default=None,
                          help="""Only generate data for these srcDbs. Can be specified multiple times.""")
        parser.add_option("--onlyDestDb", action="append", dest="onlyDestDbs", default=None,
                          help="""Only generate data for these destDbs. Can be specified multiple times.""")
        parser.add_option("--annSetType", action="append", dest="annSetTypes", default=None,
                          help="""Only generate data for these annotation sets. Can be specified multiple times, valid values are: """ + ", ".join([str(t) for t in genomeData.AnnotationType]))

    def __init__(self, opts):
        "construct GenomeDefs and mappings"
        self.defs = genomeData.load(opts.genomeDataPickle)

        # restrictions; default to full set of items to make code simpler
        self.onlySrcDbs = self.defs.mapDbNamesSet(opts.onlySrcDbs) if (opts.onlySrcDbs != None) else None
        self.onlyDestDbs = self.defs.mapDbNamesSet(opts.onlyDestDbs) if (opts.onlyDestDbs != None) else None
        self.otherDestDbs = self.defs.mapDbNamesSet(opts.otherDestDbs)
        self.excludeDbs = self.defs.mapDbNamesSet(opts.excludeDbs)
        self.clades = setOps.mkfzset(opts.includeClades) if (opts.includeClades != None) else self.defs.clades

        if opts.annotationTypes != None:
            # convert to Enumeration
            self.annotationTypes = set()
            for c in opts.annotationTypes:
                self.annotationTypes.add(genomeData.AnnotationType(c))
        else:
            self.annotationTypes = self.defs.annotationTypes
        self.annotationTypes = frozenset(self.annotationTypes)

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

    def __getDestDbsForOrgs(self):
        destDbs = set()
        for org in self.defs.orgs.itervalues():
            if self.__inclDestOrg(org):
                db = self.__getNewestDestDb(org)
                if db != None:
                    destDbs.add(db)
        return destDbs

    def __getDestDbs(self):
        "get set of destDbs, given restrictions"
        destDbs = set(self.otherDestDbs)
        if self.onlyDestDbs != None:
            destDbs |= self.onlyDestDbs
        else:
            destDbs |= self.__getDestDbsForOrgs()
        return destDbs

    def __useMapping(self, chainsSet):
        # make sure we have `all' chains (this was missing once), as we fall back
        # to using these and syntenic filtering.
        return ((chainsSet.srcDb not in self.excludeDbs) and chainsSet.haveChains()
                and (genomeData.ChainType.all in chainsSet.byType)
                and chainsSet.srcDb.matches(self.clades, self.annotationTypes))

    def __inclSrcDb(self, srcDb):
        return (self.onlySrcDbs == None) or (srcDb in self.onlySrcDbs)

    def __getSrcOrgDestDbMappings(self, srcOrg, destDb, annotationType, chainsSets, gms):
        # chainsSets are sorted new to oldest srcDb
        for chainsSet in chainsSets:
            if self.__useMapping(chainsSet) and self.__inclSrcDb(chainsSet.srcDb) and  (annotationType in chainsSet.srcDb.annotationTypes):
                gms.addMapping(chainsSet, annotationType)
                break # stop at newest

    def __getDestDbMappings(self, destDb, gms):
        "build all GenomeMapping objects for a destDb"
        for srcOrg in destDb.srcChainsSets.iterkeys():
            for annotationType in self.annotationTypes:
                self.__getSrcOrgDestDbMappings(srcOrg, destDb, annotationType, destDb.srcChainsSets[srcOrg], gms)

    def __getMappingSets(self):
        """Get a GenomeMappings object with chains describing possible
        mappings given the restrictions."""
        gms = GenomeMappings()
        for destDb in self.__getDestDbs():
            self.__getDestDbMappings(destDb, gms)
        return gms

class GenomeMapping(object):
    "mapping of a cDNA type from a source to desc db"
    def __init__(self, chainsSet, annotationType):
        self.srcDb = chainsSet.srcDb
        self.annotationType = annotationType
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

    def __writeSrcDbMappings(self, fh, srcDb):
        gms = self.mappingsBySrcDb[srcDb]
        annotationTypes = frozenset([gm.annotationType for gm in gms])
        chainTypes = frozenset(gms[0].chainsSet.byType.iterkeys())
        fileOps.prRowv(fh, srcDb.name, self.destDb.name, setOps.setJoin(annotationTypes, ","),
                       setOps.setJoin(chainTypes, ","), gms[0].chainsSet.dist)
        
    mappingsTsvHeader = ("srcDb", "destDb", "annotationTypes", "chainTypes", "dist")

    def writeMappingsTsv(self, fh):
        for srcDb in sorted(self.srcDbs, key=lambda a: a.name):
            self.__writeSrcDbMappings(fh, srcDb)
        

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
    
    def addMapping(self, chainsSet, annotationType):
        "add a mapping"
        destDbMappings = self.__obtainDestDb(chainsSet.destDb)
        gm = GenomeMapping(chainsSet, annotationType)
        destDbMappings.addMapping(gm)
        self.bySrcDb.add(gm.srcDb, gm)

    def writeMappingsTsv(self, fh):
        fileOps.prRow(fh, DestDbMappings.mappingsTsvHeader)
        for destDb in sorted(self.destDbMappingsTbl.iterkeys(), key=lambda a: a.name):
            self.destDbMappingsTbl[destDb].writeMappingsTsv(fh)
