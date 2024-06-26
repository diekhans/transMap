PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcMetadata (
            srcId text not null,
            accv text not null,
            cds text,
            geneName text,
            geneId text,
            geneType text,
            transcriptType text);
INSERT INTO srcMetadata VALUES('hg19:ENST00000254770.2','ENST00000254770.2','579..1931','LANCL2','ENSG00000132434.5','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000285279.5','ENST00000285279.5','202..720','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000414113.1','ENST00000414113.1','284..539','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000417399.1','ENST00000417399.1','434..591','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000418904.1','ENST00000418904.1','198..665','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000427700.1','ENST00000427700.1','42..554','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000428097.1','ENST00000428097.1','436..753','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000428648.1','ENST00000428648.1','308..625','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000433959.1','ENST00000433959.1','168..659','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000452107.2','ENST00000452107.2','1..122','LANCL2','ENSG00000132434.5','protein_coding','nonsense_mediated_decay');
INSERT INTO srcMetadata VALUES('hg19:ENST00000452832.1','ENST00000452832.1','421..485','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000453112.1','ENST00000453112.1','131..208','VOPP1','ENSG00000154978.8','protein_coding','nonsense_mediated_decay');
INSERT INTO srcMetadata VALUES('hg19:ENST00000453256.1','ENST00000453256.1','231..548','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000454227.1','ENST00000454227.1','135..464','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000454777.1','ENST00000454777.1',NULL,'RP11-310H4.1','ENSG00000223475.1','antisense','antisense');
INSERT INTO srcMetadata VALUES('hg19:ENST00000455023.1','ENST00000455023.1','287..574','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:ENST00000462326.1','ENST00000462326.1',NULL,'VOPP1','ENSG00000154978.8','protein_coding','processed_transcript');
INSERT INTO srcMetadata VALUES('hg19:ENST00000466041.1','ENST00000466041.1',NULL,'LANCL2','ENSG00000132434.5','protein_coding','retained_intron');
INSERT INTO srcMetadata VALUES('hg19:ENST00000471168.1','ENST00000471168.1',NULL,'VOPP1','ENSG00000154978.8','protein_coding','processed_transcript');
INSERT INTO srcMetadata VALUES('hg19:ENST00000476479.1','ENST00000476479.1',NULL,'LANCL2','ENSG00000132434.5','protein_coding','processed_transcript');
INSERT INTO srcMetadata VALUES('hg19:ENST00000481197.1','ENST00000481197.1',NULL,'VOPP1','ENSG00000154978.8','protein_coding','retained_intron');
INSERT INTO srcMetadata VALUES('hg19:ENST00000486376.1','ENST00000486376.1',NULL,'LANCL2','ENSG00000132434.5','protein_coding','processed_transcript');
INSERT INTO srcMetadata VALUES('hg19:ENST00000545390.1','ENST00000545390.1','1..510','VOPP1','ENSG00000154978.8','protein_coding','protein_coding');
COMMIT;
