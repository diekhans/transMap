PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcSeq (
            name text not null,
            seq text not null);
INSERT INTO "srcSeq" VALUES('hg19:AA004506.1','ggaattacccatcatcattaagcaatgaaacagaccggctggtgcactggtncacggcgccccgggggtcatccacatgctcatgcaggcgtacaaggtctttaaggaggagaagtacttgaaagaggccatggagtgtagcgatgtgatttggcagcgaggtttgctgcggaagggctacgggatatgccatgggactgctggcaacggctattccttcctgtccctttaccgtctcacgcaggataagaagtacctctaccgagcttgcaacgtttgcagagtggtgtctanattaccggagcacacgggtccngtattcctgacagaccctattcgctctttgaaggcatggctggcgctattcactttctctctgatgtcctgggaccagagacatcacggtttccagcatttgaacttgactcttccgaagagggattaaaaggtgcaaaaagacaactaanatacccatttggaccaa');
INSERT INTO "srcSeq" VALUES('hg19:AA004507.1','gtgaagaagtcatattttatattgatgagggtgctgttagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtnacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagccgccancaatgccttcaaagagcgaatagggtctgtcaaggaatacggcacccgtgtgctccgtaatctagacaccactctgcaaacttncaagcttcggtagaggtnacttcttatcctgcgtgagaacgggtaaag');
INSERT INTO "srcSeq" VALUES('hg19:AA021160.1','gcgctggagaggacgcgaggagccatgagngccactgcngaaggtggcggcgtgctgctcgggctgctcttggagtgcacagaagccaaaaagcattgctggtatttcgaaggactctatccaacctattatatatgccgctcctacgaggactgctgtggctccaggtgctgtgtgcgggccctctccatacagaggctgtggtacttctggttccttctgatgatgggcgtgcttttctgctgcggancggcttcttcatccggaggcgcatgtaccccccgccgctgatcgaggagccagccttcaatgtgtcctacaccaaggcagcccccaaatcccggcccaggaagcccaagcaagccggggccgccctatttacacttgaccccaggaaggaccggggatgaaccttgtcggggaattccatgggcaatgggctttt');
INSERT INTO "srcSeq" VALUES('hg19:AA318041.1','ggctgtggtacttntggttccttctgatgatgggcgtgcttttttgctgcggancggcttcttcatccggaggcgnatgtaccccccgccgctgatcgaggagccagccttcantgtgtcctacaccaggcagnccccaaatcccggcccaggagcccagcagccggggccgccctattacaccgacccaggaggaccggggatgaaccctgtcgggaattccatggcaatggctttccaggtcccacccaactcaccccaggggagt');
INSERT INTO "srcSeq" VALUES('hg19:AA488805.1','caaaagtggaccaagaaaccttgacagaaatggtgaacccagtattgattatgtgcgccacaaaaaattccgatctgggaattacccatcatcattaagcaatgaaacagaccggctggtgcactggtgccacggcgccccgggggtcatccacatgctcatgcaggcgtacaaggtctttaaggaggagaagtacttgaaagaggccatggagtgtagcgatgtgatttggcagcgaggtttgctgcggaagggctacgggatatgc');
INSERT INTO "srcSeq" VALUES('hg19:AI078217.1','ttttttttgtgaagaagtcatattttatattgatgagggtgctgttagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgcaagctcggtagaggtacttctta');
INSERT INTO "srcSeq" VALUES('hg19:AI141190.1','tctgtatgtgaagaagtcatattttatattgatgagggtgctgttagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgcaagctc');
INSERT INTO "srcSeq" VALUES('hg19:AI189888.1','tgtatgtgaagaagtcatatttgatattgatgagggtgctgttagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgcaagctcggtagaggtacttcttatcctgcgtgagacggtaaagggacaggaag');
INSERT INTO "srcSeq" VALUES('hg19:AI333456.1','tttgacaacaaacaccaagaaacatgcaaacactacagagaatctgtatgtgaagaagtcatattttatattgatgagggtgctgttagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgca');
INSERT INTO "srcSeq" VALUES('hg19:AI672671.1','gcgaagaagtcatattttatattgatgagggtgctgctagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctagtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgcaagctcggtagaggtacttcttatcctgcgtgagacggtaaagggacaggaaggaatagccgttgccagcagtcccatggcatatcccgtagcccttccgcagcaaacctcgctgccaaatcacatcgctacactccatggcctctttcaagtacttctcctccttaaagaccttgtacgcctgcatgagcatgtggatgaccccccggggcgc');
INSERT INTO "srcSeq" VALUES('hg19:AI873495.1','ttttatgtgaagaagtcatggtttatattgatgagggtgctgttagaaaatgttgatggatggcttcggtcactcattctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtttcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgcaagctcggtagaggtacttcttatcctgcgtgagacggtaaagggacaggaaggaatagccgttgccagcagtcccatggcatatcccgtagcccttccgagcanacctcgctgccaatcacatcgctacactcatggncctntcangtactttctcctcttaagaccttgtacgctgcctcgtgc');
INSERT INTO "srcSeq" VALUES('hg19:AI910316.1','tggatgcaaaagcggaccaagaaaccttgacagaaatggtgaaacccagtattgattatgtgcgccataaaaaattccgatctgggaattacccatcatcattaagcaatgaaacagaccggctggtgcactggtgccacggcgccccgggggtcatccacatgctcatgcaggcgtacaaggtctttaaggaggagaagtacttgaaagaggccatggagtgtagcgatgtgacttggcagcgaggattgctgcggaagggctacgggatatgccatgggactgctcgcaacggctattcct');
INSERT INTO "srcSeq" VALUES('hg19:AI936654.1','tttgtatgtgaagaagtcatattttatattgatgagggggggggtagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagttgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcaccttttaatccctcttcgaagagtcaagttcaaatgctggaaaccgtgatgtctttggtcccagacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactctgcaaacttgcaagctcggtagaggtacttcttatcctgcgtgagacggtaaagggacaggaaggaatagccgttgccagcagtcccatggcatatcccgtagcccttccgcagcaaacctcgctgccaaatcacatcgctacactccatggcctctttcaagtacttctcctccttanagaccttgtacgcctgcatgagcatgtggatgacccccggngcgccgtggcaccagtgcaccagcccggtctgttcattgcttaatgatgatggganattnccagatcggaattttttgtggcgcacataatcaatactggggttcaccatttctgtcaa');
INSERT INTO "srcSeq" VALUES('hg19:AI951073.1','gagttcatatattggactccatggaaagcctgaaagagagctgtgcttgctgtgaggatatcagaggaactgcccttagcagcccacgagaccgttcctggaagtgaacatcaacgaagacagaaaggccagggaaaggccctctcctgtctctcctcttgcacgtgggcaccccactacttggccttcactacctgttcgtacgggggcggaggcgtgttgcagtaggctggagggggcgggcaggccacactcccctggggtgagttgggtgggacctggaaagccattgccatggaattcccgacagggttcatccccggtcctcctgggtcggtgtaatagggcggccccggctgctgggctcctgggccgggatttgggggctgcctggtgtaggacacattgaaggctggctcctcgatcagcg');
INSERT INTO "srcSeq" VALUES('hg19:AL556266.2','gtgctgttagaaaatgttgatggatggcttcggtcactcatgctagtctaaccagaggggcctgtgacggtgtctgcttctctttcaggattcccagtcgtttctgtgtcaggcactaagcaatctggcggcttttggtccaaatgggtattttagttgtctttttgcacctcttaatccctctccgaagagtcaagtncaaatgctggaaaccgtgatgtctctggtcccaggacatcagagagaaagtgaatagcgccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgccccgtaatctagacaccactctgcaaacttgcaagcccggtagaggtacttcttatcctgcgtgagagggtaaagggacaggaaggaatagccgttgccagcagtcccatggcatatcccgtagcccttccgcagcaaacctcgctgccaaatcacatcgctacactccatggcctctttcaagtacttctcctccttaaagaccttgtacgcctgcctgagcatgtggatgacccccggggcgccgtggcaccagtncaccagccggtctgtctcattgcttaatgatgatgggtaattcccagatcggaattttttgtggcgcacataatcaatactgggtttcaccatttctgtcaaggtttcttggtccacttttgctgccggctgcattaacatatagtaaattccagccatgccatgggctgctccaacgtactgcttccggtgccactggtacaacagcgggcagcgctccgtttttctttcttcccttgacaaagtcttacccgattcaataatagcattgactacctctttaatagctgactcacacacggtgcctggacctatctctgtgttcaggtacagtnaggcatacagataacctgcccgtccataaagcagcccatcaggaaggtctgatttttggcagaaaccgatctctggagccgcaaagttttgtgacaaatnctgggactanagtacttctgagtgtgtaaaatcannccccacagcagggggcagctcccncagagg');
INSERT INTO "srcSeq" VALUES('hg19:AL556267.3','gggatgcagaaatggaggaacgggcgtncgtcaaccccttcccggactacgaggccgccgccggggcgctgctcgcctccggagcggccgaagagacaggctgtgttcgtcccccggcgaccacggatgagcccggcctcccttttcatcaggacgggaagatcattcataatttcataagacggatccagaccaaaattaaagatcttctgcagcaaatggaagaagggctgaagacagctgatccccatgactgctctgcttatactggctggacaggcatagcccttttgtacctgcagttgtaccgggtcacatgtgaccaaacctacctgctccgatccctggattacgtaaaaagaacacttcggaatctgaatggccgcaggtcaccttcctctgtggggatgctggccccctggctgttggagctgtgatttatcacaaactcagaagtgactgtgagtcccaggaatgtgtcacaaaacttttgcagctccagagatcggttgtctgccaagaatcagaccttcctgatgagctgctttatggacgggcaggttatctgtatgccttactgtacctgaacacagagataggtccaggcaccgtgtgttagtnagctattaaagaggtagtcaatgctattattgaatcgggtaagactttgtcaagggaagaaagaaaaacggagcgctgcccgctgttgtaccagtggcaccggaagcagtacgttggagcagcccatggcatggctggaatttactatatgttaatgcagccggcagcaaaagtggaccaagaaaccttgacagaaatggtgaacccagtattgattatgtgcgccacaaaaaattccgatctgggaattacccatcatcattaagcaatgaaacagaccggctggtgcactggtgccacggcgccccgggggtcatccacatgctcatgcaggcgtacaaggtctttaaggaggagaagtacttgaaagangccatggagtgtagcgatgtgatttggcagcnaggtttgctg');
INSERT INTO "srcSeq" VALUES('hg19:AL558370.3','ccgggatcccancggtcggcgggacggctcccggctgcagtctgcccgcccgccccgcgcgggggccgagtcgcgaagcgcctgcgacccggcgtccgggcgcgctggagaggacgcgaggagccatgaggcgccagctgcgaaggtggcggcgctgctgctcgggctgctcttggagtgcacagaagccaaaaagcattgctggtatttcgaaggactctatccaacctattatatatgccgctcctacgaggactgctgtggctccaggtactgtgtgcgggccctctccatacagaggctgtggtacttctggttccttctgatgatgggcgtgcttttctgctgcggancggcttcttcatccggaggcgcatgtaccccccgccgctgatcgaggagccagccttcaatgtgtcctacaccaggcagcccccaaatcccggcccaggagcccagcagccggggccgccctattacaccgacccaggaggaccggggatgaaccctgtcgggaattccatggcaatggctttccaggtcccacccaactcaccccaggggagtgtggcctgcccgccccctccagcctactgcaacacgcctccgcccccgtacgaacaggtagtgaaggccaagtagtggggtgcccacgtgcaagaggagagacaggagagggcctttccctggcctttctgtcttcgttgatgttcacttccaggaacggtctcgtgggctgctaagggcagttcctctgatatcctcacagcaagcacagctctctttcaggctttccatggagtacaatatatgaactcacactttgtctcctctgttgcttctgtttctgacgcagctggtgctctcacatggtagtgtggtgacagtccccgagggctgacgtccttacggtggcgtgaccagatctacaggagagagactgagagganaangcatnctggagtgcagtggcatgtagaggggc');
INSERT INTO "srcSeq" VALUES('hg19:AL559143.3','cgcgcgggggcgagtcgcganggcctgcgacccggcgtccgggcgcgctggagaggacgcgaggagccatgaggcgccagctgcgaaggtggcgagcgctgctgctcgggctgctcttggagtgcacagaagccaaaaagcattgctggtatttcgaaggactctatccaacctattatatatgccgctcctacgaggactgctgtggctccaggtgctgtgtgcgggccctctccatacagaggctgtggtacttctggttccttctgatgatgggcgtgcttttctgctgcggancggcttcttcatccggaggcgcatgtaccccccgccgctgatcgaggagccagccttcaatgtgtcctacaccaggcagcccccaaatcccggcccaggagcccagcagccggggccgccctattacaccgacccaggaggaccggggatgaaccctgtcgggaattccatggcaatggctttccaggnnccacccaactnaccccaggggagtgtggcctgcccgccccctccagcctactgcaacacgcctccgcccccgtacgaacaggtagtgaaggccaagtagtggggtgcccacgtgcaagaggagagacaggagagggcctttccctggcctttctgtcttcgttgatgttcacttccaggaacggtctcgtgggctgctaagggcagttcctctgatatcctcacagcaagcacagctctctttcaggctttccatggagtacaatatatganctcacactttgtctcctctgttgcttctgtttctgacgcanctgatgctctcacatggtagtgtggtgacagtccccgagggctgacgtccttacggtggcgtgaccagatctacagg');
INSERT INTO "srcSeq" VALUES('hg19:AU280098.1','agtctgcccgcccgccccgcgcgggggccgagtcgcgaagcgcgcctgcgacccggcgtccgggcgcgctggagaggacgctgaggagccatgaggcgccagcctgcgaaggtggcnggcgctgctgctcgggctgctcttggagtgcacagaagccaaaaagcattgctggtatttcgaaggactctatccaacctattatatatgccgctcctacgaggactgctgtggctccaggtgctgtgtgcgggccctctccatacagaggctgtggtacttctggttccttctgatgatgggcgtgcttttctgctgcggagccggcttcttcatccggaggcgcatgtaccccccgccgctgatcgaggagccagccttcaatgtgtcctacaccaggcagcccccaaatcccggcccaggagcccagcagccggggccgccctattacaccgacccaggaggaccggggatgaaccctgtcgggaattccatggcaatggctttccaggtcccacccaactcaccccaggggagtgtggcctg');
INSERT INTO "srcSeq" VALUES('hg19:AW261983.1','gttcaaatgctggaaaccgtgatgtctttggtcccaggacattaaagagaaagtgaatagggccagccatgccttcaaagagcgaatagggtctgtcaggaatgcggcacccgtgtgctccgtaatctagacaccactttgcaaacttgcaagctcggtggaggtactttttatcctgcttgagacggtaaaaggacaggaaggaatatccgttgccaggagtcccatggcatatgtccgaggcctt');
CREATE UNIQUE INDEX srcSeq_name on srcSeq (name);
COMMIT;
