PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE genomeAsms (
            hgDb text not null,
            clade text not null,
            commonName text not null,
            scientificName text not null,
            annotationTypeSet text);
INSERT INTO "genomeAsms" VALUES('hg19','mammal','Human','Homo sapiens','rna,est,refseq,gencode');
INSERT INTO "genomeAsms" VALUES('mm10','mammal','Mouse','Mus musculus','rna,est,refseq,gencode');
CREATE TABLE chains (
            srcHgDb text not null,
            destHgDb text not null,
            chainType text not null,
            chainFile text not null,
            netFile text not null);
INSERT INTO "chains" VALUES('mm10','hg19','all','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.all.chain.gz','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.net.gz');
INSERT INTO "chains" VALUES('mm10','hg19','rbest','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.rbest.chain.gz','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.rbest.net.gz');
INSERT INTO "chains" VALUES('mm10','hg19','syn','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.all.chain.gz','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.net.gz');
INSERT INTO "chains" VALUES('hg19','mm10','all','input/mm10.hg19.chain','input/mm10.hg19.net');
INSERT INTO "chains" VALUES('hg19','mm10','syn','input/mm10.hg19.chain','input/mm10.hg19.net');
CREATE INDEX genomeAsms_hgDb on genomeAsms (hgDb);
CREATE INDEX genomeAsms_commonName on genomeAsms (commonName);
CREATE INDEX genomeAsms_annotationTypeSet on genomeAsms (annotationTypeSet);
CREATE INDEX chains_srcHgDb on chains (srcHgDb);
CREATE INDEX chains_destHgDb on chains (destHgDb);
COMMIT;
