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
INSERT INTO srcMetadata VALUES('hg19:NM_001284282.1','NM_001284282.1','119..586','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001284283.1','NM_001284283.1','11..520','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001284284.1','NM_001284284.1','190..681','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321242.1','NM_001321242.1','64..702','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321243.1','NM_001321243.1','310..855','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321244.1','NM_001321244.1','275..787','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321246.1','NM_001321246.1','197..637','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321247.1','NM_001321247.1','197..622','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321248.1','NM_001321248.1','135..464','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321249.1','NM_001321249.1','198..704','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321250.1','NM_001321250.1','220..726','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321251.1','NM_001321251.1','206..712','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_018697.3','NM_018697.3','579..1931','LANCL2','55915','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_030796.5','NM_030796.5','179..697','VOPP1','81552','protein_coding','protein_coding');
CREATE UNIQUE INDEX srcMetadata_srcId on srcMetadata (srcId);
COMMIT;
