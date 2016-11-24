PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcMetadata (
            srcId text not null,
            accv text not null,
            cds text,
            geneId text,
            geneName text,
            geneType text,
            transcriptType text);
INSERT INTO "srcMetadata" VALUES('hg19:AB035966.1','AB035966.1','56..1408',NULL,'TASP','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AB097004.1','AB097004.1','56..574',NULL,NULL,'protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF353942.1','AF353942.1','56..1408',NULL,'LANCL2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF395824.1','AF395824.1','75..581',NULL,'GASP','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AJ278245.1','AJ278245.1','187..1539',NULL,'lancl2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AK000539.1','AK000539.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AK023882.1','AK023882.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AK075391.1','AK075391.1','131..649',NULL,NULL,'protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AK090825.1','AK090825.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AK092992.1','AK092992.1','168..659',NULL,NULL,'protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AK095857.1','AK095857.1','579..1931',NULL,NULL,'protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AK098425.1','AK098425.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AK124176.1','AK124176.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AK125579.1','AK125579.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AK126848.1','AK126848.1',NULL,NULL,NULL,'unknown','unknown');
CREATE UNIQUE INDEX srcMetadata_srcId on srcMetadata (srcId);
COMMIT;
