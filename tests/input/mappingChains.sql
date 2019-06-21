PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE mappingChains (
            bin int unsigned not null,
            qName text not null,
            qStart int unsigned not null,
            qEnd int unsigned not null,
            chainId int unsigned not null,
            offset int unsigned not null,
            length int unsigned not null);
INSERT INTO mappingChains VALUES(15,'chr7',55432374,55641209,704,0,10016);
INSERT INTO mappingChains VALUES(1009,'chr7',55654493,55661372,2698,10016,1767);
INSERT INTO mappingChains VALUES(1008,'chr7',55492943,55493130,1998953,11783,96);
INSERT INTO mappingChains VALUES(1008,'chr7',55479559,55479784,2187843,11879,117);
INSERT INTO mappingChains VALUES(1008,'chr7',55466108,55466327,2573090,11996,96);
INSERT INTO mappingChains VALUES(1008,'chr7',55492911,55493128,2629554,12092,112);
INSERT INTO mappingChains VALUES(1009,'chr7',55655435,55655547,2640962,12204,114);
INSERT INTO mappingChains VALUES(1008,'chr7',55496060,55496143,2933974,12318,95);
INSERT INTO mappingChains VALUES(1008,'chr7',55468910,55469018,2968400,12413,102);
INSERT INTO mappingChains VALUES(1008,'chr7',55459457,55459605,3057463,12515,121);
CREATE INDEX mappingChains_chrom_bin on mappingChains (qName, bin);
COMMIT;
