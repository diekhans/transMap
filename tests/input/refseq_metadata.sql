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
INSERT INTO srcMetadata VALUES('hg19:NM_001284282.2','NM_001284282.2','112..579','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001284283.2','NM_001284283.2','11..520','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001284284.1','NM_001284284.1','190..681','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321242.1','NM_001321242.1','64..702','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321243.2','NM_001321243.2','292..837','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321244.2','NM_001321244.2','257..769','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321246.2','NM_001321246.2','179..619','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321247.2','NM_001321247.2','179..604','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321248.2','NM_001321248.2','135..464','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321249.2','NM_001321249.2','198..704','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321250.2','NM_001321250.2','220..726','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_001321251.1','NM_001321251.1','206..712','VOPP1','81552','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_018697.4','NM_018697.4','690..2042','LANCL2','55915','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('hg19:NM_030796.5','NM_030796.5','179..697','VOPP1','81552','protein_coding','protein_coding');
COMMIT;
