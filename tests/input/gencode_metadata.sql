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
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000257724.3','ENST00000257724.3','264..1331','ENSG00000135272.5','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000284601.3','ENST00000284601.3','70..3438','ENSG00000154415.3','PPP1R3A','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000284602.1','ENST00000284602.1','1..225','ENSG00000154415.3','PPP1R3A','protein_coding','nonsense_mediated_decay');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000350908.4','ENST00000350908.4','317..2464','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000360232.4','ENST00000360232.4','13..1311','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000365464.1','ENST00000365464.1',NULL,'ENSG00000202334.1','5S_rRNA.194','rRNA','rRNA');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000378237.3','ENST00000378237.3','53..1150','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000390668.3','ENST00000390668.3','1..1371','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000393486.1','ENST00000393486.1','591..1331','ENSG00000135272.5','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000393489.3','ENST00000393489.3','383..2254','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000393491.3','ENST00000393491.3','1..1593','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000393494.2','ENST00000393494.2','280..2427','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000393500.3','ENST00000393500.3','821..1969','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000403559.4','ENST00000403559.4','375..2573','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000408937.3','ENST00000408937.3','375..2597','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000412402.1','ENST00000412402.1','387..650','ENSG00000128573.17','FOXP2','protein_coding','nonsense_mediated_decay');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000415146.1','ENST00000415146.1',NULL,'ENSG00000224595.1','AC073626.2','processed_transcript','processed_transcript');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000419883.1','ENST00000419883.1',NULL,'ENSG00000233607.1','AC068610.5','processed_transcript','processed_transcript');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000423503.1','ENST00000423503.1','181..312','ENSG00000135272.5','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000427207.1','ENST00000427207.1','1..471','ENSG00000135272.5','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000431629.1','ENST00000431629.1','1..471','ENSG00000135272.5','MDFIC','protein_coding','nonsense_mediated_decay');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000440349.1','ENST00000440349.1','420..668','ENSG00000128573.17','FOXP2','protein_coding','nonsense_mediated_decay');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000448022.1','ENST00000448022.1','184..315','ENSG00000135272.5','MDFIC','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000449795.1','ENST00000449795.1','477..1049','ENSG00000154415.3','PPP1R3A','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000459666.1','ENST00000459666.1',NULL,'ENSG00000128573.17','FOXP2','protein_coding','processed_transcript');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000462331.1','ENST00000462331.1','375..578','ENSG00000128573.17','FOXP2','protein_coding','protein_coding');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000495516.1','ENST00000495516.1',NULL,'ENSG00000128573.17','FOXP2','protein_coding','processed_transcript');
INSERT INTO "srcMetadata" VALUES('hg19:ENST00000498196.1','ENST00000498196.1','179..584','ENSG00000135272.5','MDFIC','protein_coding','protein_coding');
CREATE UNIQUE INDEX srcMetadata_srcId on srcMetadata (srcId);
COMMIT;
