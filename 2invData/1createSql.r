s <- r"(SET ECHO OFF NEWP 0 SPA 0 PAGES 99 FEED OFF TRIMS ON TAB OFF feedback off colsep '|' trimspool on hea on lin 9999
SPOOL TabNR.csv
select * from TabNR where rownum < 2;
SET ECHO OFF NEWP 0 SPA 0 PAGES 0 FEED OFF TRIMS ON TAB OFF feedback off colsep '|' trimspool on hea on
select * from TabNR;
SPOOL OFF
)"

#Get positions of plots is missing intentionally as I don't have access to them
#in the database
unlink("/tmp/sqlOra.sql")
for(tab in c("bam", "prf", "tfl")) {
  for(nr in c(3:7,9)) {
    cat(gsub("NR", nr, gsub("Tab", tab, s)), sep="\n", file = "/tmp/sqlOra.sql", append = TRUE)
  }
}

#sqlplus
#@sqlOra.sql

#bzip2 *.csv
