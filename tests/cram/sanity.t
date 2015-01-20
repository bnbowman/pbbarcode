  $ export INH5=`python -c "from pbcore import data ; print data.getCmpH5()"`
  $ export INBH51=`python -c "from pbcore import data ; print data.getBasH5s()[0]"`
  $ export INBH52=`python -c "from pbcore import data ; print data.getBasH5s()[1]"`
  $ export BARCODE_FASTA=$TESTDIR/../../etc/barcode.fasta
  $ echo $INBH51 > bas.fofn
  $ echo $INBH52 >> bas.fofn
  $ pbbarcode labelZmws $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --old --scoreMode paired $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --old --scoreMode paired --scoreFirst $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --old --scoreMode paired --scoreFirst --adapterSidePad 0 --insertSidePad 0 $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --scoreMode paired $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --scoreMode paired --scoreFirst $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --scoreMode paired --scoreFirst --adapterSidePad 0 --insertSidePad 0 $BARCODE_FASTA bas.fofn
  $ pbbarcode emitFastqs --fasta bas.fofn barcode.fofn
  $ pbbarcode emitFastqs --trim 20 bas.fofn barcode.fofn
  $ pbbarcode emitFastqs --subreads --trim 20 bas.fofn barcode.fofn
  $ cp $INH5 ./aligned_reads.cmp.h5         
  $ chmod 766 ./aligned_reads.cmp.h5
  $ pbbarcode labelAlignments barcode.fofn aligned_reads.cmp.h5  
Check that same holes get the same barcode (consistent scoring)
  $ cmph5tools.py stats --what "(Movie,HoleNumber,Barcode,AverageBarcodeScore)" aligned_reads.cmp.h5 | uniq
                         Movie                      Barcode                    AverageBarcodeScore                    HoleNumber
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   4.00                    3008                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   6.50                    2001                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   8.00                    4009                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   6.14                    2008                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   5.33                    3006                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   5.00                    1000                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   7.00                    4004                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                    bc5--bc10                                   5.00                    1006                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   4.33                    4006                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   4.00                    2006                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   4.00                    3002                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   4.00                    2006                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   4.83                    1009                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   4.00                    3002                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   4.83                    1009                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   4.00                    1000                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   5.33                    1007                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   6.50                    9                             
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   5.50                    1004                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   5.00                    2002                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   6.80                    2004                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   4.50                    4007                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   6.80                    2004                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   5.18                    3008                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   5.42                    2009                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   6.75                    2007                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   6.14                    2008                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   7.00                    1002                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   7.00                    1008                          
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                    bc5--bc10                                   6.50                    9                             
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                  14.00                    2000                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   5.67                    9                             
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                  14.00                    2000                          
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                     bc3--bc4                                   5.67                    9                             
  m110818_075520_42141_c100129202555500000315043109121112_s2_p0                     bc3--bc4                                   5.67                    8                             
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0                    bc5--bc10                                  10.00                    2003                          
