  $ export INH5=`python -c "from pbcore import data ; print data.getCmpH5()"`
  $ export INBH51=`python -c "from pbcore import data ; print data.getBasH5s()[0]"`
  $ export INBH52=`python -c "from pbcore import data ; print data.getBasH5s()[1]"`
  $ export BARCODE_FASTA=$TESTDIR/../../etc/barcode.fasta
  $ echo $INBH51 > bas.fofn
  $ echo $INBH52 >> bas.fofn
  $ pbbarcode labelZmws $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --scoreMode paired $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --scoreMode paired --scoreFirst $BARCODE_FASTA bas.fofn
  $ pbbarcode labelZmws --scoreMode paired --scoreFirst --adapterSidePad 0 --insertSidePad 0 $BARCODE_FASTA bas.fofn
  $ pbbarcode emitFastqs --fasta bas.fofn barcode.fofn
  $ pbbarcode emitFastqs --trim 20 bas.fofn barcode.fofn
  $ pbbarcode emitFastqs --subreads --trim 20 bas.fofn barcode.fofn
  $ cp $INH5 ./aligned_reads.cmp.h5         
  $ chmod 766 ./aligned_reads.cmp.h5
  $ pbbarcode labelAlignments barcode.fofn aligned_reads.cmp.h5  
