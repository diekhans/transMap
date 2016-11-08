PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcAligns (
            bin int unsigned not null,
            matches int unsigned not null,
            misMatches int unsigned not null,
            repMatches int unsigned not null,
            nCount int unsigned not null,
            qNumInsert int unsigned not null,
            qBaseInsert int unsigned not null,
            tNumInsert int unsigned not null,
            tBaseInsert int unsigned not null,
            strand text not null,
            qName text not null,
            qSize int unsigned not null,
            qStart int unsigned not null,
            qEnd int unsigned not null,
            tName text not null,
            tSize int unsigned not null,
            tStart int unsigned not null,
            tEnd int unsigned not null,
            blockCount int unsigned not null,
            blockSizes blob not null,
            qStarts text not null,
            tStarts text not null);
INSERT INTO "srcAligns" VALUES(585,1657,0,0,0,0,0,2,884,'+','hg19:ENST00000456328.2-1',1657,0,1657,'chr1',249250621,11868,14409,3,'359,109,1189,','0,359,468,','11868,12612,13220,');
INSERT INTO "srcAligns" VALUES(585,1653,0,0,0,0,0,2,888,'+','hg19:ENST00000515242.2-1',1653,0,1653,'chr1',249250621,11871,14412,3,'356,109,1188,','0,356,465,','11871,12612,13224,');
INSERT INTO "srcAligns" VALUES(585,1483,0,0,0,0,0,3,1053,'+','hg19:ENST00000518655.2-1',1483,0,1483,'chr1',249250621,11873,14409,4,'354,127,253,749,','0,354,481,734,','11873,12594,13402,13660,');
INSERT INTO "srcAligns" VALUES(585,712,0,0,0,0,0,2,832,'+','hg19:ENST00000473358.1-1',712,0,712,'chr1',249250621,29553,31097,3,'486,104,122,','0,486,590,','29553,30563,30975,');
INSERT INTO "srcAligns" VALUES(585,535,0,0,0,0,0,1,308,'+','hg19:ENST00000469289.1-1',535,0,535,'chr1',249250621,30266,31109,2,'401,134,','0,401,','30266,30975,');
INSERT INTO "srcAligns" VALUES(585,138,0,0,0,0,0,0,0,'+','hg19:ENST00000408384.1-1',138,0,138,'chr1',249250621,30365,30503,1,'138,','0,','30365,');
INSERT INTO "srcAligns" VALUES(585,1187,0,0,0,0,0,2,341,'-','hg19:ENST00000417324.1-1',1187,0,1187,'chr1',249250621,34553,36081,3,'621,205,361,','0,621,826,','34553,35276,35720,');
INSERT INTO "srcAligns" VALUES(585,590,0,0,0,0,0,1,239,'-','hg19:ENST00000461467.1-1',590,0,590,'chr1',249250621,35244,36073,2,'237,353,','0,237,','35244,35720,');
INSERT INTO "srcAligns" VALUES(585,883,0,0,0,0,0,1,3065,'+','hg19:ENST00000534990.1-1',883,0,883,'chr1',249250621,65881,69829,2,'90,793,','0,90,','65881,69036,');
INSERT INTO "srcAligns" VALUES(585,918,0,0,0,0,0,0,0,'+','hg19:ENST00000335137.3-1',918,0,918,'chr1',249250621,69090,70008,1,'918,','0,','69090,');
INSERT INTO "srcAligns" VALUES(585,2748,0,0,0,0,0,3,28890,'-','hg19:ENST00000466430.1-1',2748,0,2748,'chr1',249250621,89294,120932,4,'2335,150,105,158,','0,2335,2485,2590,','89294,92090,112699,120774,');
INSERT INTO "srcAligns" VALUES(585,1319,0,0,0,0,0,1,236,'-','hg19:ENST00000495576.1-1',1319,0,1319,'chr1',249250621,89550,91105,2,'500,819,','0,500,','89550,90286,');
INSERT INTO "srcAligns" VALUES(585,491,0,0,0,0,0,3,36497,'-','hg19:ENST00000477740.1-1',491,0,491,'chr1',249250621,92229,129217,4,'11,105,212,163,','0,11,116,328,','92229,112699,120720,129054,');
INSERT INTO "srcAligns" VALUES(585,629,0,0,0,0,0,2,17592,'-','hg19:ENST00000471248.1-1',629,0,629,'chr1',249250621,110952,129173,3,'405,105,119,','0,405,510,','110952,112699,129054,');
INSERT INTO "srcAligns" VALUES(73,336,0,0,0,0,0,1,4150,'-','hg19:ENST00000453576.2-1',336,0,336,'chr1',249250621,129080,133566,2,'143,193,','0,143,','129080,133373,');
INSERT INTO "srcAligns" VALUES(586,402,0,0,0,0,0,5,1240,'+','hg19:ENST00000544713.1-1',402,0,402,'chr1',249250621,135956,137598,6,'72,105,57,78,61,29,','0,72,177,234,312,373,','135956,136030,137354,137415,137497,137569,');
INSERT INTO "srcAligns" VALUES(586,323,0,0,0,0,0,1,227,'-','hg19:ENST00000493797.1-1',323,0,323,'chr1',249250621,139789,140339,2,'58,265,','0,58,','139789,140074,');
INSERT INTO "srcAligns" VALUES(586,4860,0,0,0,0,0,1,3374,'-','hg19:ENST00000484859.1-1',4860,0,4860,'chr1',249250621,141473,149707,2,'1538,3322,','0,1538,','141473,146385,');
INSERT INTO "srcAligns" VALUES(586,518,0,0,0,0,0,2,3506,'-','hg19:ENST00000490997.1-1',518,0,518,'chr1',249250621,142807,146831,3,'204,124,190,','0,204,328,','142807,146385,146641,');
INSERT INTO "srcAligns" VALUES(586,104,0,0,0,0,0,0,0,'-','hg19:ENST00000410691.1-1',104,0,104,'chr1',249250621,157783,157887,1,'104,','0,','157783,');
INSERT INTO "srcAligns" VALUES(586,457,0,0,0,0,0,1,623,'+','hg19:ENST00000496488.1-1',457,0,457,'chr1',249250621,160445,161525,2,'245,212,','0,245,','160445,161313,');
INSERT INTO "srcAligns" VALUES(586,972,0,0,0,0,0,5,8488,'-','hg19:ENST00000466557.1-1',972,0,972,'chr1',249250621,164404,173864,6,'387,59,66,216,132,112,','0,387,446,512,728,860,','164404,165883,168099,169048,172556,173752,');
INSERT INTO "srcAligns" VALUES(586,278,0,0,0,0,0,2,2601,'-','hg19:ENST00000491962.1-1',278,0,278,'chr1',249250621,165888,168767,3,'54,66,158,','0,54,120,','165888,168099,168609,');
INSERT INTO "srcAligns" VALUES(586,1292,0,0,0,0,0,1,5578,'-','hg19:ENST00000442116.1-1',1292,0,1292,'chr1',249250621,227614,234484,2,'1162,130,','0,1162,','227614,234354,');
INSERT INTO "srcAligns" VALUES(586,1925,0,0,0,0,0,2,8739,'-','hg19:ENST00000448958.1-1',1925,0,1925,'chr1',249250621,227894,238558,3,'882,902,141,','0,882,1784,','227894,237053,238417,');
INSERT INTO "srcAligns" VALUES(586,363,0,0,0,0,0,0,0,'-','hg19:ENST00000424429.1-1',363,0,363,'chr1',249250621,228291,228654,1,'363,','0,','228291,');
INSERT INTO "srcAligns" VALUES(73,2513,0,0,0,0,0,3,28885,'-','hg19:ENST00000424587.1-1',2513,0,2513,'chr1',249250621,235855,267253,4,'2100,150,105,158,','0,2100,2250,2355,','235855,238417,259016,267095,');
INSERT INTO "srcAligns" VALUES(586,510,0,0,0,0,0,1,1344,'-','hg19:ENST00000335577.4-1',510,0,510,'chr1',249250621,257267,259121,2,'405,105,','0,405,','257267,259016,');
INSERT INTO "srcAligns" VALUES(587,468,0,0,0,0,0,1,10177,'+','hg19:ENST00000426316.1-1',468,0,468,'chr1',249250621,317810,328455,2,'323,145,','0,323,','317810,328310,');
INSERT INTO "srcAligns" VALUES(587,573,0,0,0,0,0,2,3727,'+','hg19:ENST00000423728.1-1',573,0,573,'chr1',249250621,320161,324461,3,'492,58,23,','0,492,550,','320161,324287,324438,');
INSERT INTO "srcAligns" VALUES(587,575,0,0,0,0,0,2,320,'+','hg19:ENST00000432964.1-1',575,0,575,'chr1',249250621,320161,321056,3,'492,58,25,','0,492,550,','320161,320880,321031,');
INSERT INTO "srcAligns" VALUES(587,555,0,0,0,0,0,2,2152,'+','hg19:ENST00000440038.1-1',555,0,555,'chr1',249250621,322036,324743,3,'192,58,305,','0,192,250,','322036,324287,324438,');
INSERT INTO "srcAligns" VALUES(587,547,0,0,0,0,0,1,1677,'+','hg19:ENST00000419160.1-1',547,0,547,'chr1',249250621,322731,324955,2,'339,208,','0,339,','322731,324747,');
INSERT INTO "srcAligns" VALUES(587,1170,0,0,0,0,0,2,289,'+','hg19:ENST00000534867.1-1',1170,0,1170,'chr1',249250621,324437,325896,3,'249,406,515,','0,249,655,','324437,324718,325381,');
INSERT INTO "srcAligns" VALUES(587,1158,0,0,0,0,0,1,1180,'+','hg19:ENST00000456623.1-1',1158,0,1158,'chr1',249250621,324514,326852,2,'365,793,','0,365,','324514,326059,');
INSERT INTO "srcAligns" VALUES(587,2661,0,0,0,0,0,1,1037,'+','hg19:ENST00000425496.2-1',2661,0,2661,'chr1',249250621,324755,328453,2,'1759,902,','0,1759,','324755,327551,');
INSERT INTO "srcAligns" VALUES(587,402,0,0,0,0,0,5,459,'-','hg19:ENST00000534873.1-1',402,0,402,'chr1',249250621,326536,327397,6,'29,61,78,57,105,72,','0,29,90,168,225,330,','326536,326576,326641,326723,327218,327325,');
INSERT INTO "srcAligns" VALUES(587,336,0,0,0,0,0,1,4152,'+','hg19:ENST00000431812.1-1',336,0,336,'chr1',249250621,329783,334271,2,'193,143,','0,193,','329783,334128,');
INSERT INTO "srcAligns" VALUES(73,413,0,0,0,0,0,2,111614,'+','hg19:ENST00000455207.1-1',413,0,413,'chr1',249250621,334128,446155,3,'169,102,142,','0,169,271,','334128,439466,446013,');
INSERT INTO "srcAligns" VALUES(587,573,0,0,0,0,0,1,8094,'+','hg19:ENST00000455464.1-1',573,0,573,'chr1',249250621,334139,342806,2,'158,415,','0,158,','334139,342391,');
INSERT INTO "srcAligns" VALUES(587,995,0,0,0,0,0,0,0,'+','hg19:ENST00000426406.1-1',995,0,995,'chr1',249250621,367639,368634,1,'995,','0,','367639,');
INSERT INTO "srcAligns" VALUES(588,462,0,0,0,0,0,2,13896,'+','hg19:ENST00000440163.1-1',462,0,462,'chr1',249250621,439364,453722,3,'204,180,78,','0,204,384,','439364,446013,453644,');
INSERT INTO "srcAligns" VALUES(588,498,0,0,0,0,0,1,2558,'+','hg19:ENST00000453935.1-1',498,0,498,'chr1',249250621,450886,453942,2,'200,298,','0,200,','450886,453644,');
INSERT INTO "srcAligns" VALUES(588,406,0,0,0,0,0,1,326,'+','hg19:ENST00000431321.1-1',406,0,406,'chr1',249250621,453216,453948,2,'102,304,','0,102,','453216,453644,');
INSERT INTO "srcAligns" VALUES(588,607,0,0,0,0,0,1,6241,'-','hg19:ENST00000450983.1-1',607,0,607,'chr1',249250621,453632,460480,2,'534,73,','0,534,','453632,460407,');
INSERT INTO "srcAligns" VALUES(588,426,0,0,0,0,0,1,6213,'-','hg19:ENST00000412666.1-1',426,0,426,'chr1',249250621,453826,460465,2,'340,86,','0,340,','453826,460379,');
INSERT INTO "srcAligns" VALUES(588,676,0,0,0,0,0,2,1623,'+','hg19:ENST00000441866.1-1',676,0,676,'chr1',249250621,459655,461954,3,'337,135,204,','0,337,472,','459655,461153,461750,');
INSERT INTO "srcAligns" VALUES(588,842,0,0,0,0,0,2,1623,'-','hg19:ENST00000417636.1-1',842,0,842,'chr1',249250621,521368,523833,3,'370,135,337,','0,370,505,','521368,522200,523496,');
INSERT INTO "srcAligns" VALUES(73,607,0,0,0,0,0,1,6533,'+','hg19:ENST00000423796.1-1',607,0,607,'chr1',249250621,523008,530148,2,'73,534,','0,73,','523008,529614,');
INSERT INTO "srcAligns" VALUES(73,402,0,0,0,0,0,1,6505,'+','hg19:ENST00000450696.1-1',402,0,402,'chr1',249250621,523047,529954,2,'62,340,','0,62,','523047,529614,');
INSERT INTO "srcAligns" VALUES(589,437,0,0,0,0,0,1,326,'+','hg19:ENST00000440196.1-1',437,0,437,'chr1',249250621,529832,530595,2,'304,133,','0,304,','529832,530462,');
INSERT INTO "srcAligns" VALUES(589,498,0,0,0,0,0,1,2542,'+','hg19:ENST00000357876.4-1',498,0,498,'chr1',249250621,529838,532878,2,'298,200,','0,298,','529838,532678,');
INSERT INTO "srcAligns" VALUES(73,413,0,0,0,0,0,2,118352,'-','hg19:ENST00000440200.1-1',413,0,413,'chr1',249250621,536815,655580,3,'142,102,169,','0,142,244,','536815,543334,655411,');
INSERT INTO "srcAligns" VALUES(589,802,0,0,0,0,0,2,832,'-','hg19:ENST00000452176.1-1',802,0,802,'chr1',249250621,562756,564390,3,'447,263,92,','0,447,710,','562756,563340,564298,');
INSERT INTO "srcAligns" VALUES(589,79,0,0,0,0,0,0,0,'-','hg19:ENST00000459059.1-1',79,0,79,'chr1',249250621,566186,566265,1,'79,','0,','566186,');
INSERT INTO "srcAligns" VALUES(589,995,0,0,0,0,0,0,0,'-','hg19:ENST00000332831.2-1',995,0,995,'chr1',249250621,621058,622053,1,'995,','0,','621058,');
INSERT INTO "srcAligns" VALUES(73,629,0,0,0,0,0,2,17586,'-','hg19:ENST00000441245.1-1',629,0,629,'chr1',249250621,637315,655530,3,'405,105,119,','0,405,510,','637315,639064,655411,');
INSERT INTO "srcAligns" VALUES(73,274,0,0,0,0,0,1,16242,'-','hg19:ENST00000448605.1-1',274,0,274,'chr1',249250621,639064,655580,2,'105,169,','0,105,','639064,655411,');
INSERT INTO "srcAligns" VALUES(73,480,0,0,0,0,0,2,16030,'-','hg19:ENST00000419394.1-1',480,0,480,'chr1',249250621,639064,655574,3,'105,212,163,','0,105,317,','639064,647090,655411,');
INSERT INTO "srcAligns" VALUES(73,750,0,0,0,0,0,1,8109,'-','hg19:ENST00000414688.1-1',750,0,750,'chr1',249250621,646721,655580,2,'581,169,','0,581,','646721,655411,');
INSERT INTO "srcAligns" VALUES(590,336,0,0,0,0,0,1,4157,'-','hg19:ENST00000447954.1-1',336,0,336,'chr1',249250621,655437,659930,2,'143,193,','0,143,','655437,659737,');
INSERT INTO "srcAligns" VALUES(590,402,0,0,0,0,0,5,362,'+','hg19:ENST00000545502.1-1',402,0,402,'chr1',249250621,662322,663086,6,'72,105,57,78,61,29,','0,72,177,234,312,373,','662322,662396,662842,662903,662985,663057,');
INSERT INTO "srcAligns" VALUES(590,4860,0,0,0,0,0,1,3344,'-','hg19:ENST00000416385.1-1',4860,0,4860,'chr1',249250621,677192,685396,2,'1538,3322,','0,1538,','677192,682074,');
INSERT INTO "srcAligns" VALUES(590,104,0,0,0,0,0,0,0,'-','hg19:ENST00000411249.1-1',104,0,104,'chr1',249250621,693612,693716,1,'104,','0,','693612,');
INSERT INTO "srcAligns" VALUES(590,295,0,0,0,0,0,1,5599,'-','hg19:ENST00000417659.1-1',295,0,295,'chr1',249250621,694411,700305,2,'92,203,','0,92,','694411,700102,');
INSERT INTO "srcAligns" VALUES(590,456,0,0,0,0,0,1,623,'+','hg19:ENST00000422528.1-1',456,0,456,'chr1',249250621,696290,697369,2,'244,212,','0,244,','696290,697157,');
INSERT INTO "srcAligns" VALUES(590,1317,0,0,0,0,0,6,12453,'-','hg19:ENST00000428504.1-1',1317,0,1317,'chr1',249250621,700236,714006,7,'391,59,66,216,132,110,343,','0,391,450,516,732,864,974,','700236,701708,703927,704876,708355,709550,713663,');
INSERT INTO "srcAligns" VALUES(590,566,0,0,0,0,0,1,2844,'+','hg19:ENST00000457084.1-1',566,0,566,'chr1',249250621,714161,717571,2,'311,255,','0,311,','714161,717316,');
INSERT INTO "srcAligns" VALUES(590,441,0,0,0,0,0,2,25379,'+','hg19:ENST00000429505.1-1',441,0,441,'chr1',249250621,714435,740255,3,'37,304,100,','0,37,341,','714435,739298,740155,');
INSERT INTO "srcAligns" VALUES(590,513,0,0,0,0,0,1,2233,'+','hg19:ENST00000434264.1-1',513,0,513,'chr1',249250621,717324,720070,2,'192,321,','0,192,','717324,719749,');
INSERT INTO "srcAligns" VALUES(590,1194,0,0,0,0,0,0,0,'+','hg19:ENST00000358533.2-1',1194,0,1194,'chr1',249250621,721319,722513,1,'1194,','0,','721319,');
INSERT INTO "srcAligns" VALUES(590,523,0,0,0,0,0,3,8760,'-','hg19:ENST00000447500.1-1',523,0,523,'chr1',249250621,736258,745541,4,'285,93,50,95,','0,285,378,428,','736258,741178,743953,745446,');
INSERT INTO "srcAligns" VALUES(590,417,0,0,0,0,0,1,2250,'+','hg19:ENST00000443772.1-1',417,0,417,'chr1',249250621,740178,742845,2,'168,249,','0,168,','740178,742596,');
INSERT INTO "srcAligns" VALUES(590,581,0,0,0,0,0,2,4548,'+','hg19:ENST00000412115.1-1',581,0,581,'chr1',249250621,740311,745440,3,'35,107,439,','0,35,142,','740311,742596,745001,');
INSERT INTO "srcAligns" VALUES(590,402,0,0,0,0,0,1,7202,'-','hg19:ENST00000435300.1-1',402,0,402,'chr1',249250621,745488,753092,2,'62,340,','0,62,','745488,752752,');
INSERT INTO "srcAligns" VALUES(590,1944,0,0,0,0,0,1,520,'+','hg19:ENST00000326734.1-1',1944,0,1944,'chr1',249250621,752750,755214,2,'832,1112,','0,832,','752750,754102,');
INSERT INTO "srcAligns" VALUES(590,1317,0,0,0,0,0,0,0,'-','hg19:ENST00000473798.1-1',1317,0,1317,'chr1',249250621,761585,762902,1,'1317,','0,','761585,');
INSERT INTO "srcAligns" VALUES(590,1298,0,0,0,0,0,0,0,'-','hg19:ENST00000536430.1-1',1298,0,1298,'chr1',249250621,761588,762886,1,'1298,','0,','761588,');
INSERT INTO "srcAligns" VALUES(73,1571,0,0,0,0,0,4,25233,'+','hg19:ENST00000445118.1-1',1571,0,1571,'chr1',249250621,762987,789791,5,'168,102,184,96,1021,','0,168,270,454,550,','762987,764382,787306,788050,788770,');
INSERT INTO "srcAligns" VALUES(590,750,0,0,0,0,0,2,13322,'+','hg19:ENST00000441765.1-1',750,0,750,'chr1',249250621,763046,777118,3,'109,102,539,','0,109,211,','763046,764382,776579,');
INSERT INTO "srcAligns" VALUES(73,874,0,0,0,0,0,5,25080,'+','hg19:ENST00000449005.1-1',874,0,874,'chr1',249250621,763052,789006,6,'103,102,153,184,96,236,','0,103,205,358,542,638,','763052,764382,783033,787306,788050,788770,');
INSERT INTO "srcAligns" VALUES(73,612,0,0,0,0,0,4,24456,'+','hg19:ENST00000416570.1-1',612,0,612,'chr1',249250621,763078,788146,5,'77,102,153,184,96,','0,77,179,332,516,','763078,764382,783033,787306,788050,');
INSERT INTO "srcAligns" VALUES(590,894,0,0,0,0,0,1,9522,'+','hg19:ENST00000415295.1-1',894,0,894,'chr1',249250621,766984,777400,2,'73,821,','0,73,','766984,776579,');
INSERT INTO "srcAligns" VALUES(73,600,0,0,0,0,0,4,5185,'+','hg19:ENST00000448975.1-1',600,0,600,'chr1',249250621,783053,788838,5,'133,119,184,96,68,','0,133,252,436,532,','783053,784863,787306,788050,788770,');
INSERT INTO "srcAligns" VALUES(591,845,0,0,0,0,0,1,560,'+','hg19:ENST00000425657.1-1',845,0,845,'chr1',249250621,786727,788132,2,'763,82,','0,763,','786727,788050,');
INSERT INTO "srcAligns" VALUES(591,1807,0,0,0,0,0,2,7026,'-','hg19:ENST00000446136.1-1',1807,0,1807,'chr1',249250621,803450,812283,3,'605,1044,158,','0,605,1649,','803450,809491,812125,');
INSERT INTO "srcAligns" VALUES(591,504,0,0,0,0,0,3,1458,'-','hg19:ENST00000432963.1-1',504,0,504,'chr1',249250621,803619,805581,4,'291,49,48,116,','0,291,340,388,','803619,804006,804907,805465,');
INSERT INTO "srcAligns" VALUES(591,446,0,0,0,0,0,2,7954,'-','hg19:ENST00000427857.1-1',446,0,446,'chr1',249250621,803782,812182,3,'273,116,57,','0,273,389,','803782,810419,812125,');
INSERT INTO "srcAligns" VALUES(591,111,0,0,0,0,0,0,0,'-','hg19:ENST00000408219.1-1',111,0,111,'chr1',249250621,808846,808957,1,'111,','0,','808846,');
INSERT INTO "srcAligns" VALUES(591,543,0,0,0,0,0,2,80,'-','hg19:ENST00000539392.1-1',543,0,543,'chr1',249250621,844532,845155,3,'332,90,121,','0,332,422,','844532,844904,845034,');
INSERT INTO "srcAligns" VALUES(591,3043,0,0,0,0,0,1,471,'+','hg19:ENST00000448179.1-1',3043,0,3043,'chr1',249250621,846814,850328,2,'39,3004,','0,39,','846814,847324,');
INSERT INTO "srcAligns" VALUES(591,443,0,0,0,0,0,1,358,'+','hg19:ENST00000398216.2-1',443,0,443,'chr1',249250621,849550,850351,2,'274,169,','0,274,','849550,850182,');
INSERT INTO "srcAligns" VALUES(591,1389,0,0,0,0,0,3,1434,'-','hg19:ENST00000417705.1-1',1389,0,1389,'chr1',249250621,852249,855072,4,'851,89,91,358,','0,851,940,1031,','852249,853401,854204,854714,');
INSERT INTO "srcAligns" VALUES(591,630,0,0,0,0,0,1,535,'-','hg19:ENST00000432961.1-1',630,0,630,'chr1',249250621,852749,853914,2,'117,513,','0,117,','852749,853401,');
INSERT INTO "srcAligns" VALUES(591,3442,0,0,0,0,0,0,0,'-','hg19:ENST00000534924.1-1',3442,0,3442,'chr1',249250621,852954,856396,1,'3442,','0,','852954,');
INSERT INTO "srcAligns" VALUES(591,626,0,0,0,0,0,6,13786,'+','hg19:ENST00000420190.1-1',626,0,626,'chr1',249250621,860259,874671,7,'69,92,182,51,125,90,17,','0,69,161,343,394,519,609,','860259,861301,865534,866418,871151,874419,874654,');
INSERT INTO "srcAligns" VALUES(591,387,0,0,0,0,0,4,10257,'+','hg19:ENST00000437963.1-1',387,0,387,'chr1',249250621,860529,871173,5,'40,92,182,51,22,','0,40,132,314,365,','860529,861301,865534,866418,871151,');
INSERT INTO "srcAligns" VALUES(591,2551,0,0,0,0,0,13,16287,'+','hg19:ENST00000342066.3-1',2551,0,2551,'chr1',249250621,861117,879955,14,'63,92,182,51,125,90,186,163,116,79,500,125,111,668,','0,63,155,337,388,513,603,789,952,1068,1147,1647,1772,1883,','861117,861301,865534,866418,871151,874419,874654,876523,877515,877789,877938,878632,879077,879287,');
INSERT INTO "srcAligns" VALUES(591,2191,0,0,0,0,0,11,12073,'+','hg19:ENST00000341065.4-1',2191,0,2191,'chr1',249250621,865691,879955,12,'25,51,125,90,138,163,116,79,500,125,111,668,','0,25,76,201,291,429,592,708,787,1287,1412,1523,','865691,866418,871151,874419,874654,876523,877515,877789,877938,878632,879077,879287,');
INSERT INTO "srcAligns" VALUES(591,1731,0,0,0,0,0,6,3254,'+','hg19:ENST00000455979.1-1',1731,0,1731,'chr1',249250621,874654,879639,7,'186,163,116,79,500,125,562,','0,186,349,465,544,1044,1169,','874654,876523,877515,877789,877938,878632,879077,');
CREATE TABLE srcXRefs (
            srcAlnId text not null,
            srcId text not null,
            accv text);
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000456328.2-1','hg19:ENST00000456328.2','ENST00000456328.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000515242.2-1','hg19:ENST00000515242.2','ENST00000515242.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000518655.2-1','hg19:ENST00000518655.2','ENST00000518655.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000473358.1-1','hg19:ENST00000473358.1','ENST00000473358.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000469289.1-1','hg19:ENST00000469289.1','ENST00000469289.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000408384.1-1','hg19:ENST00000408384.1','ENST00000408384.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000417324.1-1','hg19:ENST00000417324.1','ENST00000417324.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000461467.1-1','hg19:ENST00000461467.1','ENST00000461467.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000534990.1-1','hg19:ENST00000534990.1','ENST00000534990.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000335137.3-1','hg19:ENST00000335137.3','ENST00000335137.3');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000466430.1-1','hg19:ENST00000466430.1','ENST00000466430.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000495576.1-1','hg19:ENST00000495576.1','ENST00000495576.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000477740.1-1','hg19:ENST00000477740.1','ENST00000477740.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000471248.1-1','hg19:ENST00000471248.1','ENST00000471248.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000453576.2-1','hg19:ENST00000453576.2','ENST00000453576.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000544713.1-1','hg19:ENST00000544713.1','ENST00000544713.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000493797.1-1','hg19:ENST00000493797.1','ENST00000493797.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000484859.1-1','hg19:ENST00000484859.1','ENST00000484859.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000490997.1-1','hg19:ENST00000490997.1','ENST00000490997.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000410691.1-1','hg19:ENST00000410691.1','ENST00000410691.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000496488.1-1','hg19:ENST00000496488.1','ENST00000496488.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000466557.1-1','hg19:ENST00000466557.1','ENST00000466557.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000491962.1-1','hg19:ENST00000491962.1','ENST00000491962.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000442116.1-1','hg19:ENST00000442116.1','ENST00000442116.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000448958.1-1','hg19:ENST00000448958.1','ENST00000448958.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000424429.1-1','hg19:ENST00000424429.1','ENST00000424429.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000424587.1-1','hg19:ENST00000424587.1','ENST00000424587.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000335577.4-1','hg19:ENST00000335577.4','ENST00000335577.4');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000426316.1-1','hg19:ENST00000426316.1','ENST00000426316.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000423728.1-1','hg19:ENST00000423728.1','ENST00000423728.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000432964.1-1','hg19:ENST00000432964.1','ENST00000432964.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000440038.1-1','hg19:ENST00000440038.1','ENST00000440038.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000419160.1-1','hg19:ENST00000419160.1','ENST00000419160.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000534867.1-1','hg19:ENST00000534867.1','ENST00000534867.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000456623.1-1','hg19:ENST00000456623.1','ENST00000456623.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000425496.2-1','hg19:ENST00000425496.2','ENST00000425496.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000534873.1-1','hg19:ENST00000534873.1','ENST00000534873.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000431812.1-1','hg19:ENST00000431812.1','ENST00000431812.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000455207.1-1','hg19:ENST00000455207.1','ENST00000455207.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000455464.1-1','hg19:ENST00000455464.1','ENST00000455464.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000426406.1-1','hg19:ENST00000426406.1','ENST00000426406.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000440163.1-1','hg19:ENST00000440163.1','ENST00000440163.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000453935.1-1','hg19:ENST00000453935.1','ENST00000453935.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000431321.1-1','hg19:ENST00000431321.1','ENST00000431321.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000450983.1-1','hg19:ENST00000450983.1','ENST00000450983.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000412666.1-1','hg19:ENST00000412666.1','ENST00000412666.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000441866.1-1','hg19:ENST00000441866.1','ENST00000441866.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000417636.1-1','hg19:ENST00000417636.1','ENST00000417636.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000423796.1-1','hg19:ENST00000423796.1','ENST00000423796.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000450696.1-1','hg19:ENST00000450696.1','ENST00000450696.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000440196.1-1','hg19:ENST00000440196.1','ENST00000440196.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000357876.4-1','hg19:ENST00000357876.4','ENST00000357876.4');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000440200.1-1','hg19:ENST00000440200.1','ENST00000440200.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000452176.1-1','hg19:ENST00000452176.1','ENST00000452176.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000459059.1-1','hg19:ENST00000459059.1','ENST00000459059.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000332831.2-1','hg19:ENST00000332831.2','ENST00000332831.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000441245.1-1','hg19:ENST00000441245.1','ENST00000441245.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000448605.1-1','hg19:ENST00000448605.1','ENST00000448605.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000419394.1-1','hg19:ENST00000419394.1','ENST00000419394.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000414688.1-1','hg19:ENST00000414688.1','ENST00000414688.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000447954.1-1','hg19:ENST00000447954.1','ENST00000447954.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000545502.1-1','hg19:ENST00000545502.1','ENST00000545502.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000416385.1-1','hg19:ENST00000416385.1','ENST00000416385.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000411249.1-1','hg19:ENST00000411249.1','ENST00000411249.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000417659.1-1','hg19:ENST00000417659.1','ENST00000417659.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000422528.1-1','hg19:ENST00000422528.1','ENST00000422528.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000428504.1-1','hg19:ENST00000428504.1','ENST00000428504.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000457084.1-1','hg19:ENST00000457084.1','ENST00000457084.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000429505.1-1','hg19:ENST00000429505.1','ENST00000429505.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000434264.1-1','hg19:ENST00000434264.1','ENST00000434264.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000358533.2-1','hg19:ENST00000358533.2','ENST00000358533.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000447500.1-1','hg19:ENST00000447500.1','ENST00000447500.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000443772.1-1','hg19:ENST00000443772.1','ENST00000443772.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000412115.1-1','hg19:ENST00000412115.1','ENST00000412115.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000435300.1-1','hg19:ENST00000435300.1','ENST00000435300.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000326734.1-1','hg19:ENST00000326734.1','ENST00000326734.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000473798.1-1','hg19:ENST00000473798.1','ENST00000473798.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000536430.1-1','hg19:ENST00000536430.1','ENST00000536430.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000445118.1-1','hg19:ENST00000445118.1','ENST00000445118.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000441765.1-1','hg19:ENST00000441765.1','ENST00000441765.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000449005.1-1','hg19:ENST00000449005.1','ENST00000449005.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000416570.1-1','hg19:ENST00000416570.1','ENST00000416570.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000415295.1-1','hg19:ENST00000415295.1','ENST00000415295.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000448975.1-1','hg19:ENST00000448975.1','ENST00000448975.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000425657.1-1','hg19:ENST00000425657.1','ENST00000425657.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000446136.1-1','hg19:ENST00000446136.1','ENST00000446136.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000432963.1-1','hg19:ENST00000432963.1','ENST00000432963.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000427857.1-1','hg19:ENST00000427857.1','ENST00000427857.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000408219.1-1','hg19:ENST00000408219.1','ENST00000408219.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000539392.1-1','hg19:ENST00000539392.1','ENST00000539392.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000448179.1-1','hg19:ENST00000448179.1','ENST00000448179.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000398216.2-1','hg19:ENST00000398216.2','ENST00000398216.2');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000417705.1-1','hg19:ENST00000417705.1','ENST00000417705.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000432961.1-1','hg19:ENST00000432961.1','ENST00000432961.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000534924.1-1','hg19:ENST00000534924.1','ENST00000534924.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000420190.1-1','hg19:ENST00000420190.1','ENST00000420190.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000437963.1-1','hg19:ENST00000437963.1','ENST00000437963.1');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000342066.3-1','hg19:ENST00000342066.3','ENST00000342066.3');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000341065.4-1','hg19:ENST00000341065.4','ENST00000341065.4');
INSERT INTO "srcXRefs" VALUES('hg19:ENST00000455979.1-1','hg19:ENST00000455979.1','ENST00000455979.1');
CREATE INDEX srcAligns_tName_bin on srcAligns (tName, bin);
CREATE INDEX srcAligns_qname on srcAligns (qName);
CREATE UNIQUE INDEX srcXRefs_srcAlnId on srcXRefs (srcAlnId);
CREATE INDEX srcXRefs_accv on srcXRefs (accv);
COMMIT;
