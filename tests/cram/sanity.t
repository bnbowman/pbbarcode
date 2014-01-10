  $ export INH5=`python -c "from pbcore import data ; print data.getCmpH5()['cmph5']"`
  $ export INBH51=`python -c "from pbcore import data ; print data.getCmpH5()['bash5s'][0]"`
  $ export INBH52=`python -c "from pbcore import data ; print data.getCmpH5()['bash5s'][1]"`
  $ export BARCODE_FASTA=$TESTDIR/../../etc/barcode.fasta
  $ echo $INBH51 > bas.fofn
  $ echo $INBH52 >> bas.fofn
  $ PBbarcode.py labelZmws $BARCODE_FASTA bas.fofn
  $ PBbarcode.py labelZmws --scoreMode paired $BARCODE_FASTA bas.fofn
  $ PBbarcode.py labelZmws --scoreMode paired --scoreFirst $BARCODE_FASTA bas.fofn
  $ PBbarcode.py labelZmws --scoreMode paired --scoreFirst --adapterSidePad 0 --insertSidePad 0 $BARCODE_FASTA bas.fofn
  $ PBbarcode.py emitFastqs --fasta bas.fofn barcode.fofn
  $ PBbarcode.py emitFastqs --trim 20 bas.fofn barcode.fofn
  $ PBbarcode.py emitFastqs --subreads --trim 20 bas.fofn barcode.fofn
  $ cp $INH5 ./aligned_reads.cmp.h5         
  $ chmod 766 ./aligned_reads.cmp.h5
  $ PBbarcode.py labelAlignments barcode.fofn aligned_reads.cmp.h5  
