#!/usr/bin/bash

FASTQ_DIR=
OUT_FASTQ_DIR=corrected_reads

FASTQ_HEADER_ONLY="SRR4113427_1.fastq SRR4113429_1.fastq SRR4113431_1.fastq SRR4113432_1.fastq SRR4113433_1.fastq SRR4113435_1.fastq SRR4113436_1.fastq SRR4113437_1.fastq"
FASTQ_SWAP="SRR4113428_1.fastq SRR4113430_1.fastq SRR4113434_1.fastq SRR4113438_1.fastq SRR4113439_1.fastq SRR4113440_1.fastq"

for FASTQ in $FASTQ_SWAP; do
   echo "******* Processing: $FASTQ *******"
   READS1=$FASTQ
   READS2=${FASTQ/_1.fastq/_2.fastq}

   #TPP needs header to be (slightly) different. Use the commands below to get an allowed format
   sed '1~2s/$/\/1/' $FASTQ_DIR/$READS2 > $OUT_FASTQ_DIR/$READS1
   sed '1~2s/$/\/2/' $FASTQ_DIR/$READS1 > $OUT_FASTQ_DIR/$READS2 
done

for FASTQ in $FASTQ_HEADER_ONLY; do
   echo "******* Processing: $FASTQ *******"
   READS1=$FASTQ
   READS2=${FASTQ/_1.fastq/_2.fastq}

   #TPP needs header to be (slightly) different. Use the commands below to get an allowed format
   sed '1~2s/$/\/1/' $FASTQ_DIR/$READS1 > $OUT_FASTQ_DIR/$READS1
   sed '1~2s/$/\/2/' $FASTQ_DIR/$READS2 > $OUT_FASTQ_DIR/$READS2 
done
