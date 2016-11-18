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
INSERT INTO "srcMetadata" VALUES('hg19:NM_001166345.1','NM_001166345.1','591..1331','29969','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_001166346.1','NM_001166346.1','264..722','29969','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_001172766.2','NM_001172766.2','375..2519','93986','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_001172767.2','NM_001172767.2','11..1384','93986','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_002711.3','NM_002711.3','32..3400','5506','PPP1R3A','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_014491.3','NM_014491.3','375..2522','93986','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_148898.3','NM_148898.3','375..2597','93986','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_148899.3','NM_148899.3','11..1309','93986','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_148900.3','NM_148900.3','375..2573','93986','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NM_199072.4','NM_199072.4','264..1331','29969','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NR_033766.1','NR_033766.1',NULL,'93986','FOXP2','non_coding','non_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NR_033767.1','NR_033767.1',NULL,'93986','FOXP2','non_coding','non_coding');
INSERT INTO "srcMetadata" VALUES('hg19:NR_037439.1','NR_037439.1',NULL,'100500896','MIR3666','non_coding','non_coding');
CREATE UNIQUE INDEX srcMetadata_srcId on srcMetadata (srcId);
COMMIT;
