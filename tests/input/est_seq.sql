PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE srcSeq (
            name text not null,
            seq text not null);
CREATE UNIQUE INDEX srcSeq_name on srcSeq (name);
COMMIT;
