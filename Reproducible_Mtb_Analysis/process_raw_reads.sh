#!/bin/sh

TPP=tpp

PYTHON="/usr/bin/env python2"

PATH_TO_DATA=Reads_from_SRA/corrected_reads
BWA=/usr/bin/bwa
#BWA_ALG="mem"
BWA_ALG="aln"

REFS=H37Rv_GenBank.fa
#REPLICON_ID="AL123456.3"
REPLICON_ID="AL123456"
FASTQ_DIR=$PATH_TO_DATA
OUT_DIR=intermed_output
PREFIXES_OUTFILE=$OUT_DIR/`basename $FASTQ_DIR`_prefixes.txt

# These are used for creating a CSV file
GENBANK_FILE=H37Rv_GenBank.gb
CSV_OUTFILE=H37Rv_raw_counts.csv
UNIQUE_FIELDS="locus_tag"
FIELDS="product"

PRIMER=AACCTGTTA
MISMATCHES=1
WINDOW_SIZE=6

mkdir -p $OUT_DIR

COUNTER=0
INITIAL_START_TIME=$SECONDS
#for FASTQ in $FASTQ_DIR/*_1.fastq; do
#   (( COUNTER += 1 ))
#   echo "******** Run $COUNTER: $FASTQ ********"
#   READS1=$FASTQ
#   READS2=${FASTQ/_1.fastq/_2.fastq}
#
#   OUTNAME=$(basename $FASTQ)
#   OUTNAME=${OUTNAME/_1.fastq/}
#   ITERATION_START_TIME=$SECONDS
#   $TPP -himar1 -bwa $BWA -bwa-alg $BWA_ALG -ref $REFS -replicon-ids $REPLICON_ID -reads1 $READS1 -reads2 $READS2 -window-size $WINDOW_SIZE -primer $PRIMER -mismatches $MISMATCHES -output $OUT_DIR/$OUTNAME
#   ITERATION_END_TIME=$SECONDS
#   (( ITERATION_TIME = ITERATION_END_TIME - ITERATION_START_TIME ))
#
#   (( TOTAL_RUN_TIME = SECONDS - INITIAL_START_TIME )) 
#   (( CURRENT_AVG = TOTAL_RUN_TIME / COUNTER ))
#   echo "******** TPP finished in $ITERATION_TIME seconds! Average iteration time over $COUNTER iterations:  $CURRENT_AVG seconds. ********"
#done

echo "TPP iterations done! Performing post-processing operations:"
echo "Checking TA sites hit for all samples..."
NUM_REPLICONS=`echo "$REPLICON_ID" | wc -w`
TA_HITS_OUTFILE=$OUT_DIR/`basename $FASTQ_DIR`_TAs_hit.txt
$PYTHON get_runstats.py $OUT_DIR $NUM_REPLICONS TAs_hit > $TA_HITS_OUTFILE
$PYTHON get_runstats.py $OUT_DIR $NUM_REPLICONS TAs_hit
echo ""
echo "Creating prefixes file with all prefixes from all runs..."
echo -n "" > $PREFIXES_OUTFILE
for FASTQ in $FASTQ_DIR/*_1.fastq; do
  OUTNAME=$(basename $FASTQ)
  OUTNAME=${OUTNAME/_1.fastq/}
echo $OUTNAME >> $PREFIXES_OUTFILE
done
#basename -a $OUT_DIR/*.wig | cut -c-16 | uniq > $PREFIXES_OUTFILE
echo "Created '$PREFIXES_OUTFILE'."
echo ""
echo "Creating CSV file with all samples processed by TPP..."
$PYTHON wig_gb_to_csv.py -l $PREFIXES_OUTFILE -g $GENBANK_FILE -u $UNIQUE_FIELDS -f $FIELDS -o $CSV_OUTFILE
exit
echo "Created '$CSV_OUTFILE'."
echo ""
echo "********** TPP driver script finished in a total of $TOTAL_RUN_TIME seconds **********"
