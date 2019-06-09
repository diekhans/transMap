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
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000000852.3','ENSFCAT00000000852.3','375..2852','IFT88','ENSFCAG00000000852.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000000853.2','ENSFCAT00000000853.2','1..507','EEF1AKMT1','ENSFCAG00000000853.2','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000000855.3','ENSFCAT00000000855.3','1..3399','XPO4','ENSFCAG00000000855.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000001427.3','ENSFCAT00000001427.3','447..2255',NULL,'ENSFCAG00000001427.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000002362.3','ENSFCAT00000002362.3','1..2580','MPHOSPH8','ENSFCAG00000002362.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000004915.3','ENSFCAT00000004915.3','1..4542','RNF17','ENSFCAG00000004914.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000004919.3','ENSFCAT00000004919.3','183..4202','CENPJ','ENSFCAG00000004916.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000005992.3','ENSFCAT00000005992.3','44..1225','MICU2','ENSFCAG00000005990.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000006897.3','ENSFCAT00000006897.3','1..3135','ATP12A','ENSFCAG00000006891.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000006961.2','ENSFCAT00000006961.2','1..642',NULL,'ENSFCAG00000006959.2','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000006963.2','ENSFCAT00000006963.2','1..756','GJB2','ENSFCAG00000006961.2','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000007028.2','ENSFCAT00000007028.2','231..749','SAP18','ENSFCAG00000007026.2','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000007301.3','ENSFCAT00000007301.3','102..1061','CRYL1','ENSFCAG00000007299.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000009377.3','ENSFCAT00000009377.3','1..1068','ZDHHC20','ENSFCAG00000009375.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000009836.3','ENSFCAT00000009836.3','1..1083','TNFRSF19','ENSFCAG00000009834.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000009837.3','ENSFCAT00000009837.3','964..2547','MIPEP','ENSFCAG00000009835.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000009838.1','ENSFCAT00000009838.1','1..999',NULL,'ENSFCAG00000009836.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000009839.3','ENSFCAT00000009839.3',NULL,NULL,'ENSFCAG00000009837.3','pseudogene','pseudogene');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000011122.2','ENSFCAT00000011122.2','405..3572','LATS2','ENSFCAG00000011119.2','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000021261.1','ENSFCAT00000021261.1',NULL,NULL,'ENSFCAG00000021274.1','miRNA','miRNA');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000023886.1','ENSFCAT00000023886.1','1..342',NULL,'ENSFCAG00000023124.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000023931.1','ENSFCAT00000023931.1','1..3702','PARP4','ENSFCAG00000031561.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000024557.1','ENSFCAT00000024557.1','11..1297','ZMYM5','ENSFCAG00000027258.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000024641.1','ENSFCAT00000024641.1','40..2586','RNF17','ENSFCAG00000004914.3','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000025201.1','ENSFCAT00000025201.1','2..547',NULL,'ENSFCAG00000022757.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000025509.1','ENSFCAT00000025509.1','1..1374',NULL,'ENSFCAG00000029850.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000025821.1','ENSFCAT00000025821.1','694..1320','FGF9','ENSFCAG00000030874.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000026072.1','ENSFCAT00000026072.1','1..13731','SACS','ENSFCAG00000028276.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000026124.1','ENSFCAT00000026124.1',NULL,'5S_rRNA','ENSFCAG00000023750.1','rRNA','rRNA');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000026194.1','ENSFCAT00000026194.1','1..525',NULL,'ENSFCAG00000031388.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000026389.1','ENSFCAT00000026389.1','1..1269','SKA3','ENSFCAG00000027742.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000027076.1','ENSFCAT00000027076.1','1..306','MRPL57','ENSFCAG00000024896.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000027176.1','ENSFCAT00000027176.1','1..903',NULL,'ENSFCAG00000027025.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000028000.1','ENSFCAT00000028000.1','1..627','GJB6','ENSFCAG00000021967.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000028120.1','ENSFCAT00000028120.1','1..3951','ZMYM2','ENSFCAG00000030408.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000029696.1','ENSFCAT00000029696.1','1..1038',NULL,'ENSFCAG00000021931.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000031359.1','ENSFCAT00000031359.1','1..633',NULL,'ENSFCAG00000031108.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000032059.1','ENSFCAT00000032059.1','8..883','SGCG','ENSFCAG00000029063.1','protein_coding','protein_coding');
INSERT INTO srcMetadata VALUES('felCat5:ENSFCAT00000032098.1','ENSFCAT00000032098.1','1..1566','PSPC1','ENSFCAG00000027841.1','protein_coding','protein_coding');
CREATE UNIQUE INDEX srcMetadata_srcId on srcMetadata (srcId);
COMMIT;
