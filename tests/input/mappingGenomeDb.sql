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
            srcDb text not null,
            destDb text not null,
            chainType text not null,
            chainFile text not null,
            netFile text not null,
            dist float);
INSERT INTO "chains" VALUES('mm10','hg19','all','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.all.chain.gz','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.net.gz',0.495284);
INSERT INTO "chains" VALUES('mm10','hg19','rbest','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.rbest.chain.gz','/hive/data/genomes/hg19/bed/lastz.mm10/axtChain/hg19.mm10.rbest.net.gz',0.495284);
INSERT INTO "chains" VALUES('hg19','mm10','all','input/mm10.hg19.chain','input/mm10.hg19.net',0.495284);
CREATE INDEX genomeAsms_hgDb on genomeAsms (hgDb);
CREATE INDEX genomeAsms_commonName on genomeAsms (commonName);
CREATE INDEX genomeAsms_annotationTypeSet on genomeAsms (annotationTypeSet);
CREATE INDEX chains_srcDb on chains (srcDb);
CREATE INDEX chains_destDb on chains (destDb);
COMMIT;
