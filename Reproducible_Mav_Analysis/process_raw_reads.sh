#!/bin/sh

TPP=tpp

PATH_TO_DATA=Reads_from_SRA
BWA=/usr/bin/bwa
#BWA_ALG="mem"
BWA_ALG="aln"

REFS=Avium109.fa
REPLICON_ID="CP029332,CP029333,CP029334"
FASTQ_DIR=$PATH_TO_DATA
TPP_OUT_DIR=intermed_output
PREFIXES_OUTFILE=$TPP_OUT_DIR/`basename $FASTQ_DIR`_prefixes.txt

# These are used for creating a CSV file
GENBANK_FILE=Avium109.gb
CSV_OUTFILE=MAC109_raw_counts.csv
UNIQUE_FIELDS="locus_tag"
FIELDS="product regulatory_class bound_moiety"

PRIMER=AACCTGTTA
MISMATCHES=2
WINDOW_SIZE=6

mkdir -p $TPP_OUT_DIR

#COUNTER=0
#INITIAL_START_TIME=$SECONDS
#for FASTQ in $FASTQ_DIR/*1_1.fastq; do
#  (( COUNTER += 1 ))
#  echo "******** Run $COUNTER: $FASTQ ********"
#  READS1=$FASTQ
#  READS2=${FASTQ/_1.fastq/_2.fastq}
#
#  OUTNAME=$(basename $FASTQ)
#  OUTNAME=${OUTNAME/_1.fastq/}
#  ITERATION_START_TIME=$SECONDS
#  $TPP -himar1 -bwa $BWA -bwa-alg $BWA_ALG -ref $REFS -replicon-ids $REPLICON_ID -reads1 $READS1 -reads2 $READS2 -window-size $WINDOW_SIZE -primer $PRIMER -mismatches $MISMATCHES -output $TPP_OUT_DIR/$OUTNAME
#  ITERATION_END_TIME=$SECONDS
#  (( ITERATION_TIME = ITERATION_END_TIME - ITERATION_START_TIME ))
# 
#  (( TOTAL_RUN_TIME = SECONDS - INITIAL_START_TIME )) 
#  (( CURRENT_AVG = TOTAL_RUN_TIME / COUNTER ))
#  echo "******** TPP finished in $ITERATION_TIME seconds! Average iteration time over $COUNTER iterations:  $CURRENT_AVG seconds. ********"
#done

echo "TPP iterations done! Performing post-processing operations:"
echo "Checking TA sites hit for all samples..."
NUM_REPLICONS=`echo "$REPLICON_ID" | wc -w`
TA_HITS_OUTFILE=$TPP_OUT_DIR/`basename $FASTQ_DIR`_TAs_hit.txt
python2 get_runstats.py $TPP_OUT_DIR $NUM_REPLICONS TAs_hit > $TA_HITS_OUTFILE
python2 get_runstats.py $TPP_OUT_DIR $NUM_REPLICONS TAs_hit
echo ""
echo "Creating prefixes file with all prefixes from all runs..."
echo -n "" > $PREFIXES_OUTFILE
for FASTQ in $FASTQ_DIR/*_1.fastq; do
  OUTNAME=$(basename $FASTQ)
  OUTNAME=${OUTNAME/_1.fastq/}
echo $OUTNAME >> $PREFIXES_OUTFILE
done
#basename -a $TPP_OUT_DIR/*_1.fastq ${OUTNAME/_1.fastq/} > $PREFIXES_OUTFILE
#ls -1 $TPP_OUT_DIR/*.wig | cut -c-18 | uniq > $PREFIXES_OUTFILE
echo "Created '$PREFIXES_OUTFILE'."
echo ""
echo "Creating CSV file with all samples processed by TPP..."
python2 wig_gb_to_csv.py -l $PREFIXES_OUTFILE -g $GENBANK_FILE -u $UNIQUE_FIELDS -f $FIELDS -o $CSV_OUTFILE
echo "Created '$CSV_OUTFILE'."
echo ""
(( TOTAL_RUN_TIME = SECONDS - INITIAL_START_TIME ))
echo "********** process_raw_reads.sh finished in a total of $TOTAL_RUN_TIME seconds **********"
