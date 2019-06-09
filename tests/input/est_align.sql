PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcAlign (
            bin INT UNSIGNED NOT NULL,
            matches INT UNSIGNED NOT NULL,
            misMatches INT UNSIGNED NOT NULL,
            repMatches INT UNSIGNED NOT NULL,
            nCount INT UNSIGNED NOT NULL,
            qNumInsert INT UNSIGNED NOT NULL,
            qBaseInsert INT UNSIGNED NOT NULL,
            tNumInsert INT UNSIGNED NOT NULL,
            tBaseInsert INT UNSIGNED NOT NULL,
            strand TEXT NOT NULL,
            qName TEXT NOT NULL,
            qSize INT UNSIGNED NOT NULL,
            qStart INT UNSIGNED NOT NULL,
            qEnd INT UNSIGNED NOT NULL,
            tName TEXT NOT NULL,
            tSize INT UNSIGNED NOT NULL,
            tStart INT UNSIGNED NOT NULL,
            tEnd INT UNSIGNED NOT NULL,
            blockCount INT UNSIGNED NOT NULL,
            blockSizes TEXT NOT NULL,
            qStarts TEXT NOT NULL,
            tStarts TEXT NOT NULL);
INSERT INTO srcAlign VALUES(125,1029,2,0,4,0,0,8,58222,'+','hg19:AL556267.3-1',1035,0,1035,'chr7',159138663,55433761,55493018,9,'161,118,110,97,148,147,44,138,72,','0,161,279,389,486,634,781,825,963,','55433761,55459485,55466115,55466226,55467649,55468866,55479599,55479644,55492946,');
INSERT INTO srcAlign VALUES(1008,977,12,0,3,0,0,6,30542,'-','hg19:AL556266.2-1',1073,0,992,'chr7',159138663,55467652,55499186,7,'22,122,147,183,177,73,268,','81,103,225,372,555,732,805,','55467652,55467675,55468866,55479599,55492946,55496067,55498918,');
INSERT INTO srcAlign VALUES(1008,293,5,0,0,0,0,1,13164,'+','hg19:AI910316.1-1',303,5,303,'chr7',159138663,55479605,55493067,2,'177,121,','5,182,','55479605,55492946,');
INSERT INTO srcAlign VALUES(1008,268,0,0,0,0,0,2,13165,'+','hg19:AA488805.1-1',268,0,268,'chr7',159138663,55479606,55493039,3,'37,138,93,','0,37,175,','55479606,55479644,55492946,');
INSERT INTO srcAlign VALUES(1008,696,9,0,4,1,1,5,18888,'-','hg19:AI936654.1-1',710,0,710,'chr7',159138663,55479626,55499223,7,'87,6,62,177,73,35,269,','0,87,94,156,333,406,441,','55479626,55479714,55479720,55492946,55496067,55498918,55498954,');
INSERT INTO srcAlign VALUES(1008,475,2,0,3,4,4,4,18888,'+','hg19:AA004506.1-1',484,0,484,'chr7',159138663,55479684,55499052,7,'51,45,176,24,50,82,52,','0,52,97,274,299,349,432,','55479684,55479737,55492946,55496066,55496090,55498918,55499000,');
INSERT INTO srcAlign VALUES(1008,585,4,0,0,1,1,3,18886,'-','hg19:AI672671.1-1',590,0,590,'chr7',159138663,55479741,55499216,5,'12,29,177,73,298,','0,13,42,219,292,','55479741,55479753,55492946,55496067,55498918,');
INSERT INTO srcAlign VALUES(1008,538,5,0,4,1,1,8,5728,'-','hg19:AI873495.1-1',568,3,551,'chr7',159138663,55492944,55499219,10,'5,4,8,11,6,18,16,105,73,301,','17,22,26,35,46,52,70,86,191,264,','55492944,55492950,55492955,55492963,55492975,55492982,55493001,55493018,55496067,55498918,');
INSERT INTO srcAlign VALUES(1008,224,11,0,0,0,0,2,5722,'-','hg19:AW261983.1-1',247,0,235,'chr7',159138663,55493033,55498990,3,'90,73,72,','12,102,175,','55493033,55496067,55498918,');
INSERT INTO srcAlign VALUES(1008,432,2,0,0,0,0,2,5722,'-','hg19:AI189888.1-1',434,0,434,'chr7',159138663,55493065,55499221,3,'58,73,303,','0,58,131,','55493065,55496067,55498918,');
INSERT INTO srcAlign VALUES(1008,399,1,0,0,0,0,2,5722,'-','hg19:AI078217.1-1',407,7,407,'chr7',159138663,55493095,55499217,3,'28,73,299,','0,28,101,','55493095,55496067,55498918,');
INSERT INTO srcAlign VALUES(1008,368,3,0,2,2,2,1,2778,'-','hg19:AA004507.1-1',425,0,375,'chr7',159138663,55496065,55499216,4,'49,26,11,287,','50,100,126,138,','55496065,55496114,55498918,55498929,');
INSERT INTO srcAlign VALUES(1008,379,1,0,0,0,0,1,2778,'-','hg19:AI141190.1-1',389,0,380,'chr7',159138663,55496065,55499223,2,'75,305,','9,84,','55496065,55498918,');
INSERT INTO srcAlign VALUES(1008,418,1,0,0,0,0,1,2778,'-','hg19:AI333456.1-1',426,3,422,'chr7',159138663,55496065,55499262,2,'75,344,','4,79,','55496065,55498918,');
INSERT INTO srcAlign VALUES(15,904,4,60,2,2,8,10,98967,'-','hg19:AL558370.3-1',985,7,985,'chr7',159138663,55540200,55640137,11,'17,5,7,498,98,37,78,59,41,60,70,','0,17,22,36,534,633,670,748,807,848,908,','55540200,55540218,55540224,55540240,55559974,55560074,55565305,55588764,55639963,55640005,55640067,');
INSERT INTO srcAlign VALUES(15,853,4,12,6,2,2,8,98958,'-','hg19:AL559143.3-1',877,0,877,'chr7',159138663,55540256,55640089,10,'482,98,37,78,59,27,13,59,12,10,','0,482,581,618,696,755,783,796,855,867,','55540256,55559974,55560074,55565305,55588764,55639963,55639990,55640004,55640066,55640079,');
INSERT INTO srcAlign VALUES(1008,429,1,0,0,0,0,1,19236,'+','hg19:AI951073.1-1',430,0,430,'chr7',159138663,55540371,55560037,2,'367,63,','0,367,','55540371,55559974,');
INSERT INTO srcAlign VALUES(15,511,0,30,0,2,2,4,98951,'-','hg19:AU280098.1-1',543,0,543,'chr7',159138663,55540614,55640106,7,'124,137,78,59,28,34,81,','0,124,261,339,398,427,462,','55540614,55559974,55565305,55588764,55639963,55639991,55640025,');
INSERT INTO srcAlign VALUES(1008,248,1,0,3,1,1,2,19238,'-','hg19:AA318041.1-1',268,15,268,'chr7',159138663,55540622,55560112,3,'116,98,38,','0,116,215,','55540622,55559974,55560074,');
INSERT INTO srcAlign VALUES(15,342,0,0,0,4,8,7,79726,'-','hg19:AA021160.1-1',445,0,350,'chr7',159138663,55559973,55640041,9,'26,73,37,78,59,23,13,4,29,','95,122,196,233,311,370,393,411,416,','55559973,55559999,55560074,55565305,55588764,55639963,55639987,55640006,55640012,');
CREATE TABLE srcXRef (
            srcAlignId text not null,
            srcId text not null,
            accv text);
INSERT INTO srcXRef VALUES('hg19:AL556267.3-1','hg19:AL556267.3','AL556267.3');
INSERT INTO srcXRef VALUES('hg19:AL556266.2-1','hg19:AL556266.2','AL556266.2');
INSERT INTO srcXRef VALUES('hg19:AI910316.1-1','hg19:AI910316.1','AI910316.1');
INSERT INTO srcXRef VALUES('hg19:AA488805.1-1','hg19:AA488805.1','AA488805.1');
INSERT INTO srcXRef VALUES('hg19:AI936654.1-1','hg19:AI936654.1','AI936654.1');
INSERT INTO srcXRef VALUES('hg19:AA004506.1-1','hg19:AA004506.1','AA004506.1');
INSERT INTO srcXRef VALUES('hg19:AI672671.1-1','hg19:AI672671.1','AI672671.1');
INSERT INTO srcXRef VALUES('hg19:AI873495.1-1','hg19:AI873495.1','AI873495.1');
INSERT INTO srcXRef VALUES('hg19:AW261983.1-1','hg19:AW261983.1','AW261983.1');
INSERT INTO srcXRef VALUES('hg19:AI189888.1-1','hg19:AI189888.1','AI189888.1');
INSERT INTO srcXRef VALUES('hg19:AI078217.1-1','hg19:AI078217.1','AI078217.1');
INSERT INTO srcXRef VALUES('hg19:AA004507.1-1','hg19:AA004507.1','AA004507.1');
INSERT INTO srcXRef VALUES('hg19:AI141190.1-1','hg19:AI141190.1','AI141190.1');
INSERT INTO srcXRef VALUES('hg19:AI333456.1-1','hg19:AI333456.1','AI333456.1');
INSERT INTO srcXRef VALUES('hg19:AL558370.3-1','hg19:AL558370.3','AL558370.3');
INSERT INTO srcXRef VALUES('hg19:AL559143.3-1','hg19:AL559143.3','AL559143.3');
INSERT INTO srcXRef VALUES('hg19:AI951073.1-1','hg19:AI951073.1','AI951073.1');
INSERT INTO srcXRef VALUES('hg19:AU280098.1-1','hg19:AU280098.1','AU280098.1');
INSERT INTO srcXRef VALUES('hg19:AA318041.1-1','hg19:AA318041.1','AA318041.1');
INSERT INTO srcXRef VALUES('hg19:AA021160.1-1','hg19:AA021160.1','AA021160.1');
CREATE INDEX srcAlign_tName_bin ON srcAlign (tName, bin);
CREATE INDEX srcAlign_qname ON srcAlign (qName);
CREATE UNIQUE INDEX srcXRef_srcAlignId on srcXRef (srcAlignId);
CREATE INDEX srcXRef_accv on srcXRef (accv);
COMMIT;
