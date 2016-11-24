PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE mappingChains (
            bin int unsigned not null,
            qName text not null,
            qStart int unsigned not null,
            qEnd int unsigned not null,
            offset int unsigned not null,
            length int unsigned not null);
INSERT INTO "mappingChains" VALUES(15,'chr7',55432374,55641209,0,10016);
INSERT INTO "mappingChains" VALUES(1009,'chr7',55649919,55650243,10016,129);
INSERT INTO "mappingChains" VALUES(1008,'chr7',55492911,55493128,10145,112);
CREATE INDEX mappingChains_chrom_bin on mappingChains (qName, bin);
COMMIT;
