PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcAlign (
            bin int unsigned not null,
            matches int unsigned not null,
            misMatches int unsigned not null,
            repMatches int unsigned not null,
            nCount int unsigned not null,
            qNumInsert int unsigned not null,
            qBaseInsert int unsigned not null,
            tNumInsert int unsigned not null,
            tBaseInsert int unsigned not null,
            strand text not null,
            qName text not null,
            qSize int unsigned not null,
            qStart int unsigned not null,
            qEnd int unsigned not null,
            tName text not null,
            tSize int unsigned not null,
            tStart int unsigned not null,
            tEnd int unsigned not null,
            blockCount int unsigned not null,
            blockSizes blob not null,
            qStarts text not null,
            tStarts text not null);
CREATE TABLE srcXRef (
            srcAlignId text not null,
            srcId text not null,
            accv text);
CREATE INDEX srcAlign_tName_bin on srcAlign (tName, bin);
CREATE INDEX srcAlign_qname on srcAlign (qName);
CREATE UNIQUE INDEX srcXRef_srcAlignId on srcXRef (srcAlignId);
CREATE INDEX srcXRef_accv on srcXRef (accv);
COMMIT;
