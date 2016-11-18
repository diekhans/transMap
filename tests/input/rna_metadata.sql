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
INSERT INTO "srcMetadata" VALUES('hg19:AF024579.1','AF024579.1','1..225',NULL,'PPP1R3','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF054589.1','AF054589.1','264..1331',NULL,NULL,'protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF086040.1','AF086040.1',NULL,NULL,NULL,'unknown','unknown');
INSERT INTO "srcMetadata" VALUES('hg19:AF337817.1','AF337817.1','151..2298',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF454830.1','AF454830.1','1063..>1509',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467252.1','AF467252.1','102..>766',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467253.1','AF467253.1','13..>1257',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467254.1','AF467254.1','<1..290',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467255.1','AF467255.1','<1..526',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467256.1','AF467256.1','<1..290',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467257.1','AF467257.1','53..1150',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467258.1','AF467258.1','1..>1293',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF467259.1','AF467259.1','1..>1368',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AF493430.1','AF493430.1','168..674',NULL,'FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:AK131266.1','AK131266.1','797..1909',NULL,NULL,'protein_coding','protein_coding');
CREATE UNIQUE INDEX srcMetadata_srcId on srcMetadata (srcId);
COMMIT;
