#!/bin/bash

: ${1?"Usage: $0 <output directory for mtDNA reads> <output directory for hg38 reads> <results subdirectory for source library output>"}

# Set up script
shopt -s expand_aliases
set -e

mtDNA_DIR=$1
hg38_DIR=$2
RESULTS_DIR=$3

echo "-- Generate combined mtDNA fastq files --\n"

for bam in ${RESULTS_DIR}/*/prealignments/*_rCRSd.bam
do
  # Extract sample name and path
  export NAME=$(echo "$bam" | sed -e 's|${RESULTS_DIR}/||' -e 's|/prealignments.*$||')
  export SAMPLE_PATH=$(echo "$bam" | sed -e 's|/prealignments/.*$||')
  echo "$NAME mtDNA fastq files"
  # Then convert to fastq
  if [ ! -f "${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r1.fq" ]; then
    samtools bam2fq -@ 16 $bam > ${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r1.fq
  fi
  # Re-pair
  if [ ! -f "${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r2.fq" ]; then
    cp ${SAMPLE_PATH}/raw/${NAME}_R2.fastq.gz ${SAMPLE_PATH}/prealignments/
    pigz -df ${SAMPLE_PATH}/prealignments/${NAME}_R2.fastq.gz
    fastq_pair ${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r1.fq ${SAMPLE_PATH}/prealignments/${NAME}_R2.fastq
    mv -f ${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r1.fq.paired.fq ${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r1.fq
    mv -f ${SAMPLE_PATH}/prealignments/${NAME}_R2.fastq.paired.fq ${SAMPLE_PATH}/prealignments/${NAME}_mtDNA_r2.fq
  fi
done

# Combine fastq files into total pool
if [ ! -f "${mtDNA_DIR}pooled_mtDNA_r1.fq" ]; then
    cat ${RESULTS_DIR}/*/prealignments/*_mtDNA_r1.fq > ${mtDNA_DIR}pooled_mtDNA_r1.fq
fi
if [ ! -f "${mtDNA_DIR}pooled_mtDNA_r2.fq" ]; then
    cat ${RESULTS_DIR}/*/prealignments/*_mtDNA_r2.fq > ${mtDNA_DIR}pooled_mtDNA_r2.fq
fi

echo "-- Generate combined hg38 fastq files --\n"

for bam in ${RESULTS_DIR}/*/aligned_hg38/*_sort_dedup.bam
do
  # Extract sample name and path
  export NAME=$(echo "$bam" | sed -e 's|${RESULTS_DIR}||' -e 's|/aligned_hg38.*$||')
  export SAMPLE_PATH=$(echo "$bam" | sed -e 's|/aligned_hg38/.*$||')
  echo "$NAME hg38 fastq files"
  # Then convert to fastq
  if [ ! -f "${SAMPLE_PATH}/aligned_hg38/${NAME}_hg38_r1.fq" ]; then
    samtools bam2fq -@ 16 -1 ${SAMPLE_PATH}/aligned_hg38/${NAME}_hg38_r1.fq -2 ${SAMPLE_PATH}/aligned_hg38/${NAME}_hg38_r2.fq -0 /dev/null -s /dev/null $bam
  fi
done

if [ ! -f "${hg38_DIR}pooled_hg38_r1.fq" ]; then
    cat ${RESULTS_DIR}/*/aligned_hg38/*_hg38_r1.fq > ${hg38_DIR}pooled_hg38_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_r2.fq" ]; then
    cat ${RESULTS_DIR}/*/aligned_hg38/*_hg38_r2.fq > ${hg38_DIR}pooled_hg38_r2.fq
fi

echo "-- Creating pooled read files --"

if [ ! -f "${hg38_DIR}pooled_hg38_1M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 1000000 > ${hg38_DIR}pooled_hg38_1M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_2M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 2000000 > ${hg38_DIR}pooled_hg38_2M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_3M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 3000000 > ${hg38_DIR}pooled_hg38_3M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_4M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 4000000 > ${hg38_DIR}pooled_hg38_4M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_5M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 5000000 > ${hg38_DIR}pooled_hg38_5M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_6M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 6000000 > ${hg38_DIR}pooled_hg38_6M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_7M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 7000000 > ${hg38_DIR}pooled_hg38_7M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_8M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 8000000 > ${hg38_DIR}pooled_hg38_8M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_9M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 9000000 > ${hg38_DIR}pooled_hg38_9M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_10M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 10000000 > ${hg38_DIR}pooled_hg38_10M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_1M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 1000000 > ${hg38_DIR}pooled_hg38_1M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_2M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 2000000 > ${hg38_DIR}pooled_hg38_2M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_3M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 3000000 > ${hg38_DIR}pooled_hg38_3M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_4M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 4000000 > ${hg38_DIR}pooled_hg38_4M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_5M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 5000000 > ${hg38_DIR}pooled_hg38_5M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_6M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 6000000 > ${hg38_DIR}pooled_hg38_6M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_7M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 7000000 > ${hg38_DIR}pooled_hg38_7M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_8M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 8000000 > ${hg38_DIR}pooled_hg38_8M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_9M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 9000000 > ${hg38_DIR}pooled_hg38_9M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_10M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 10000000 > ${hg38_DIR}pooled_hg38_10M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_12M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 12000000 > ${mtDNA_DIR}pooled_mtDNA_12M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_14M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 14000000 > ${mtDNA_DIR}pooled_mtDNA_14M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_16M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 16000000 > ${mtDNA_DIR}pooled_mtDNA_16M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_18M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 18000000 > ${mtDNA_DIR}pooled_mtDNA_18M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_20M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 20000000 > ${mtDNA_DIR}pooled_mtDNA_20M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_12M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 12000000 > ${mtDNA_DIR}pooled_mtDNA_12M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_14M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 14000000 > ${mtDNA_DIR}pooled_mtDNA_14M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_16M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 16000000 > ${mtDNA_DIR}pooled_mtDNA_16M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_18M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 18000000 > ${mtDNA_DIR}pooled_mtDNA_18M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_20M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 20000000 > ${mtDNA_DIR}pooled_mtDNA_20M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_12M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 12000000 > ${hg38_DIR}pooled_hg38_12M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_14M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 14000000 > ${hg38_DIR}pooled_hg38_14M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_16M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 16000000 > ${hg38_DIR}pooled_hg38_16M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_18M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 18000000 > ${hg38_DIR}pooled_hg38_18M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_20M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 20000000 > ${hg38_DIR}pooled_hg38_20M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_12M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 12000000 > ${hg38_DIR}pooled_hg38_12M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_14M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 14000000 > ${hg38_DIR}pooled_hg38_14M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_16M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 16000000 > ${hg38_DIR}pooled_hg38_16M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_18M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 18000000 > ${hg38_DIR}pooled_hg38_18M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_20M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 20000000 > ${hg38_DIR}pooled_hg38_20M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_24M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 24000000 > ${mtDNA_DIR}pooled_mtDNA_24M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_28M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 28000000 > ${mtDNA_DIR}pooled_mtDNA_28M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_32M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 32000000 > ${mtDNA_DIR}pooled_mtDNA_32M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_36M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 36000000 > ${mtDNA_DIR}pooled_mtDNA_36M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_40M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 40000000 > ${mtDNA_DIR}pooled_mtDNA_40M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_24M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 24000000 > ${mtDNA_DIR}pooled_mtDNA_24M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_28M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 28000000 > ${mtDNA_DIR}pooled_mtDNA_28M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_32M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 32000000 > ${mtDNA_DIR}pooled_mtDNA_32M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_36M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 36000000 > ${mtDNA_DIR}pooled_mtDNA_36M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_40M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 40000000 > ${mtDNA_DIR}pooled_mtDNA_40M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_24M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 24000000 > ${hg38_DIR}pooled_hg38_24M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_28M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 28000000 > ${hg38_DIR}pooled_hg38_28M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_32M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 32000000 > ${hg38_DIR}pooled_hg38_32M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_36M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 36000000 > ${hg38_DIR}pooled_hg38_36M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_40M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 40000000 > ${hg38_DIR}pooled_hg38_40M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_24M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 24000000 > ${hg38_DIR}pooled_hg38_24M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_28M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 28000000 > ${hg38_DIR}pooled_hg38_28M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_32M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 32000000 > ${hg38_DIR}pooled_hg38_32M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_36M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 36000000 > ${hg38_DIR}pooled_hg38_36M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_hg38_40M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 40000000 > ${hg38_DIR}pooled_hg38_40M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_30M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 30000000 > ${mtDNA_DIR}pooled_mtDNA_30M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_42M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 42000000 > ${mtDNA_DIR}pooled_mtDNA_42M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_48M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 48000000 > ${mtDNA_DIR}pooled_mtDNA_48M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_54M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 54000000 > ${mtDNA_DIR}pooled_mtDNA_54M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_60M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 60000000 > ${mtDNA_DIR}pooled_mtDNA_60M_r1.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_30M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 30000000 > ${mtDNA_DIR}pooled_mtDNA_30M_r2.fq
fi
if [ ! -f "${hg38_DIR}pooled_mtDNA_42M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 42000000 > ${mtDNA_DIR}pooled_mtDNA_42M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_48M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 48000000 > ${mtDNA_DIR}pooled_mtDNA_48M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_54M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 54000000 > ${mtDNA_DIR}pooled_mtDNA_54M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_60M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 60000000 > ${mtDNA_DIR}pooled_mtDNA_60M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_30M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 30000000 > ${hg38_DIR}pooled_hg38_30M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_42M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 42000000 > ${hg38_DIR}pooled_hg38_42M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_48M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 48000000 > ${hg38_DIR}pooled_hg38_48M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_54M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 54000000 > ${hg38_DIR}pooled_hg38_54M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_60M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 60000000 > ${hg38_DIR}pooled_hg38_60M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_30M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 30000000 > ${hg38_DIR}pooled_hg38_30M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_42M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 42000000 > ${hg38_DIR}pooled_hg38_42M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_48M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 48000000 > ${hg38_DIR}pooled_hg38_48M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_54M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 54000000 > ${hg38_DIR}pooled_hg38_54M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_60M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 60000000 > ${hg38_DIR}pooled_hg38_60M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_56M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 56000000 > ${mtDNA_DIR}pooled_mtDNA_56M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_64M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 64000000 > ${mtDNA_DIR}pooled_mtDNA_64M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_72M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 72000000 > ${mtDNA_DIR}pooled_mtDNA_72M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_80M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 80000000 > ${mtDNA_DIR}pooled_mtDNA_80M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_56M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 56000000 > ${mtDNA_DIR}pooled_mtDNA_56M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_64M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 64000000 > ${mtDNA_DIR}pooled_mtDNA_64M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_72M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 72000000 > ${mtDNA_DIR}pooled_mtDNA_72M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_80M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 80000000 > ${mtDNA_DIR}pooled_mtDNA_80M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_56M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 56000000 > ${hg38_DIR}pooled_hg38_56M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_64M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 64000000 > ${hg38_DIR}pooled_hg38_64M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_72M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 72000000 > ${hg38_DIR}pooled_hg38_72M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_80M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 80000000 > ${hg38_DIR}pooled_hg38_80M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_56M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 56000000 > ${hg38_DIR}pooled_hg38_56M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_64M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 64000000 > ${hg38_DIR}pooled_hg38_64M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_72M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 72000000 > ${hg38_DIR}pooled_hg38_72M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_80M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 80000000 > ${hg38_DIR}pooled_hg38_80M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_50M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 50000000 > ${mtDNA_DIR}pooled_mtDNA_50M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_70M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 70000000 > ${mtDNA_DIR}pooled_mtDNA_70M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_90M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 90000000 > ${mtDNA_DIR}pooled_mtDNA_90M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_100M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 100000000 > ${mtDNA_DIR}pooled_mtDNA_100M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_50M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 50000000 > ${mtDNA_DIR}pooled_mtDNA_50M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_70M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 70000000 > ${mtDNA_DIR}pooled_mtDNA_70M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_90M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 90000000 > ${mtDNA_DIR}pooled_mtDNA_90M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_100M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 100000000 > ${mtDNA_DIR}pooled_mtDNA_100M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_50M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 50000000 > ${hg38_DIR}pooled_hg38_50M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_70M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 70000000 > ${hg38_DIR}pooled_hg38_70M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_90M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 90000000 > ${hg38_DIR}pooled_hg38_90M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_100M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 100000000 > ${hg38_DIR}pooled_hg38_100M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_50M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 50000000 > ${hg38_DIR}pooled_hg38_50M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_70M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 70000000 > ${hg38_DIR}pooled_hg38_70M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_90M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 90000000 > ${hg38_DIR}pooled_hg38_90M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_100M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 100000000 > ${hg38_DIR}pooled_hg38_100M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_84M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 84000000 > ${mtDNA_DIR}pooled_mtDNA_84M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_96M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 96000000 > ${mtDNA_DIR}pooled_mtDNA_96M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_108M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 108000000 > ${mtDNA_DIR}pooled_mtDNA_108M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_120M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 120000000 > ${mtDNA_DIR}pooled_mtDNA_120M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_84M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 84000000 > ${mtDNA_DIR}pooled_mtDNA_84M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_96M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 96000000 > ${mtDNA_DIR}pooled_mtDNA_96M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_108M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 108000000 > ${mtDNA_DIR}pooled_mtDNA_108M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_120M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 120000000 > ${mtDNA_DIR}pooled_mtDNA_120M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_84M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 84000000 > ${hg38_DIR}pooled_hg38_84M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_96M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 96000000 > ${hg38_DIR}pooled_hg38_96M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_108M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 108000000 > ${hg38_DIR}pooled_hg38_108M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_120M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 120000000 > ${hg38_DIR}pooled_hg38_120M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_84M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 84000000 > ${hg38_DIR}pooled_hg38_84M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_96M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 96000000 > ${hg38_DIR}pooled_hg38_96M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_108M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 108000000 > ${hg38_DIR}pooled_hg38_108M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_120M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 120000000 > ${hg38_DIR}pooled_hg38_120M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_98M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 98000000 > ${mtDNA_DIR}pooled_mtDNA_98M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_112M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 112000000 > ${mtDNA_DIR}pooled_mtDNA_112M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_126M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 126000000 > ${mtDNA_DIR}pooled_mtDNA_126M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_140M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 140000000 > ${mtDNA_DIR}pooled_mtDNA_140M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_98M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 98000000 > ${mtDNA_DIR}pooled_mtDNA_98M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_112M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 112000000 > ${mtDNA_DIR}pooled_mtDNA_112M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_126M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 126000000 > ${mtDNA_DIR}pooled_mtDNA_126M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_140M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 140000000 > ${mtDNA_DIR}pooled_mtDNA_140M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_98M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 98000000 > ${hg38_DIR}pooled_hg38_98M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_112M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 112000000 > ${hg38_DIR}pooled_hg38_112M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_126M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 126000000 > ${hg38_DIR}pooled_hg38_126M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_140M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 140000000 > ${hg38_DIR}pooled_hg38_140M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_98M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 98000000 > ${hg38_DIR}pooled_hg38_98M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_112M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 112000000 > ${hg38_DIR}pooled_hg38_112M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_126M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 126000000 > ${hg38_DIR}pooled_hg38_126M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_140M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 140000000 > ${hg38_DIR}pooled_hg38_140M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_128M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 128000000 > ${mtDNA_DIR}pooled_mtDNA_128M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_144M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 144000000 > ${mtDNA_DIR}pooled_mtDNA_144M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_160M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 160000000 > ${mtDNA_DIR}pooled_mtDNA_160M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_128M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 128000000 > ${mtDNA_DIR}pooled_mtDNA_128M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_144M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 144000000 > ${mtDNA_DIR}pooled_mtDNA_144M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_160M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 160000000 > ${mtDNA_DIR}pooled_mtDNA_160M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_128M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 128000000 > ${hg38_DIR}pooled_hg38_128M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_144M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 144000000 > ${hg38_DIR}pooled_hg38_144M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_160M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 160000000 > ${hg38_DIR}pooled_hg38_160M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_128M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 128000000 > ${hg38_DIR}pooled_hg38_128M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_144M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 144000000 > ${hg38_DIR}pooled_hg38_144M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_160M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 160000000 > ${hg38_DIR}pooled_hg38_160M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_162M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 162000000 > ${mtDNA_DIR}pooled_mtDNA_162M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_180M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 180000000 > ${mtDNA_DIR}pooled_mtDNA_180M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_162M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 162000000 > ${mtDNA_DIR}pooled_mtDNA_162M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_180M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 180000000 > ${mtDNA_DIR}pooled_mtDNA_180M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_162M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 162000000 > ${hg38_DIR}pooled_hg38_162M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_180M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 180000000 > ${hg38_DIR}pooled_hg38_180M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_162M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 162000000 > ${hg38_DIR}pooled_hg38_162M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_180M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 180000000 > ${hg38_DIR}pooled_hg38_180M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_200M_r1.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r1.fq 200000000 > ${mtDNA_DIR}pooled_mtDNA_200M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_mtDNA_200M_r2.fq" ]; then
    seqtk sample -s99 ${mtDNA_DIR}pooled_mtDNA_r2.fq 200000000 > ${mtDNA_DIR}pooled_mtDNA_200M_r2.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_200M_r1.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r1.fq 200000000 > ${hg38_DIR}pooled_hg38_200M_r1.fq
fi

if [ ! -f "${hg38_DIR}pooled_hg38_200M_r2.fq" ]; then
    seqtk sample -s99 ${hg38_DIR}pooled_hg38_r2.fq 200000000 > ${hg38_DIR}pooled_hg38_200M_r2.fq
fi



function percent_total() {
    echo "-- Generate 0% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}0-10_${TOTAL}_r2.fq" ]; then
      cat ${hg38_DIR}pooled_hg38_${TOTAL}_r1.fq > ${mtDNA_DIR}0-10_${TOTAL}_r1.fq
      cat ${hg38_DIR}pooled_hg38_${TOTAL}_r2.fq > ${mtDNA_DIR}0-10_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 10% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}1-9_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt10}_r1.fq ${hg38_DIR}pooled_hg38_${hg10}_r1.fq > ${mtDNA_DIR}1-9_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt10}_r2.fq ${hg38_DIR}pooled_hg38_${hg10}_r2.fq > ${mtDNA_DIR}1-9_${TOTAL}_r2.fq
    fi
    
    echo "-- Generat20% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}2-8_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt20}_r1.fq ${hg38_DIR}pooled_hg38_${hg20}_r1.fq > ${mtDNA_DIR}2-8_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt20}_r2.fq ${hg38_DIR}pooled_hg38_${hg20}_r2.fq > ${mtDNA_DIR}2-8_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 30% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}3-7_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt30}_r1.fq ${hg38_DIR}pooled_hg38_${hg30}_r1.fq > ${mtDNA_DIR}3-7_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt30}_r2.fq ${hg38_DIR}pooled_hg38_${hg30}_r2.fq > ${mtDNA_DIR}3-7_${TOTAL}_r2.fq
    fi
    
    
    echo "-- Generate 40% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}4-6_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt40}_r1.fq ${hg38_DIR}pooled_hg38_${hg40}_r1.fq > ${mtDNA_DIR}4-6_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt40}_r2.fq ${hg38_DIR}pooled_hg38_${hg40}_r2.fq > ${mtDNA_DIR}4-6_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 50% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}5-5_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt50}_r1.fq ${hg38_DIR}pooled_hg38_${hg50}_r1.fq > ${mtDNA_DIR}5-5_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt50}_r2.fq ${hg38_DIR}pooled_hg38_${hg50}_r2.fq > ${mtDNA_DIR}5-5_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 60% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}6-4_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt60}_r1.fq ${hg38_DIR}pooled_hg38_${hg60}_r1.fq > ${mtDNA_DIR}6-4_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt60}_r2.fq ${hg38_DIR}pooled_hg38_${hg60}_r2.fq > ${mtDNA_DIR}6-4_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 70% mtDNA --"

    if [ ! -f "${mtDNA_DIR}7-3_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt70}_r1.fq ${hg38_DIR}pooled_hg38_${hg70}_r1.fq > ${mtDNA_DIR}7-3_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt70}_r2.fq ${hg38_DIR}pooled_hg38_${hg70}_r2.fq > ${mtDNA_DIR}7-3_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 80% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}8-2_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt80}_r1.fq ${hg38_DIR}pooled_hg38_${hg80}_r1.fq > ${mtDNA_DIR}8-2_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt80}_r2.fq ${hg38_DIR}pooled_hg38_${hg80}_r2.fq > ${mtDNA_DIR}8-2_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 90% mtDNA --"
    
    if [ ! -f "${mtDNA_DIR}9-1_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt90}_r1.fq ${hg38_DIR}pooled_hg38_${hg90}_r1.fq > ${mtDNA_DIR}9-1_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt90}_r2.fq ${hg38_DIR}pooled_hg38_${hg90}_r2.fq > ${mtDNA_DIR}9-1_${TOTAL}_r2.fq
    fi
    
    echo "-- Generate 100% mtDNA --"
    
    if [ ! -f "$${mtDNA_DIR}10-0_${TOTAL}_r2.fq" ]; then
      cat ${mtDNA_DIR}pooled_mtDNA_${mt100}_r1.fq > ${mtDNA_DIR}10-0_${TOTAL}_r1.fq
      cat ${mtDNA_DIR}pooled_mtDNA_${mt100}_r2.fq > ${mtDNA_DIR}10-0_${TOTAL}_r2.fq
    fi
}

echo "-- Produce 10M total pooled read sets --\n"

export TOTAL=10M

export mt10=1M
export mt20=2M
export mt30=3M
export mt40=4M
export mt50=5M
export mt60=6M
export mt70=7M
export mt80=8M
export mt90=9M
export mt100=10M
export hg10=9M
export hg20=8M
export hg30=7M
export hg40=6M
export hg50=5M
export hg60=4M
export hg70=3M
export hg80=2M
export hg90=1M

percent_total

echo "-- Produce 20M total pooled read sets --\n"

export TOTAL=20M

export mt10=2M
export mt20=4M
export mt30=6M
export mt40=8M
export mt50=10M
export mt60=12M
export mt70=14M
export mt80=16M
export mt90=18M
export mt100=20M
export hg10=18M
export hg20=16M
export hg30=14M
export hg40=12M
export hg50=10M
export hg60=8M
export hg70=6M
export hg80=4M
export hg90=2M

percent_total

echo "-- Produce 40M total pooled read sets --\n"

export TOTAL=40M

export mt10=4M
export mt20=8M
export mt30=12M
export mt40=16M
export mt50=20M
export mt60=24M
export mt70=28M
export mt80=32M
export mt90=36M
export mt100=40M
export hg10=36M
export hg20=32M
export hg30=28M
export hg40=24M
export hg50=20M
export hg60=16M
export hg70=12M
export hg80=8M
export hg90=4M

percent_total

echo "-- Produce 60M total pooled read sets --\n"

export TOTAL=60M

export mt10=6M
export mt20=12M
export mt30=18M
export mt40=24M
export mt50=30M
export mt60=36M
export mt70=42M
export mt80=48M
export mt90=54M
export mt100=60M
export hg10=54M
export hg20=48M
export hg30=42M
export hg40=36M
export hg50=30M
export hg60=24M
export hg70=18M
export hg80=12M
export hg90=6M

percent_total

echo "-- Produce 80M total pooled read sets --\n"

export TOTAL=80M

export mt10=8M
export mt20=16M
export mt30=24M
export mt40=32M
export mt50=40M
export mt60=48M
export mt70=56M
export mt80=64M
export mt90=72M
export mt100=80M
export hg10=72M
export hg20=64M
export hg30=56M
export hg40=48M
export hg50=40M
export hg60=32M
export hg70=24M
export hg80=16M
export hg90=8M

percent_total

echo "-- Produce 100M total pooled read sets --\n"

export TOTAL=100M

export mt10=10M
export mt20=20M
export mt30=30M
export mt40=40M
export mt50=50M
export mt60=60M
export mt70=70M
export mt80=80M
export mt90=90M
export mt100=100M
export hg10=90M
export hg20=80M
export hg30=70M
export hg40=60M
export hg50=50M
export hg60=40M
export hg70=30M
export hg80=20M
export hg90=10M

percent_total

echo "-- Produce 120M total pooled read sets --\n"

export TOTAL=120M

export mt10=12M
export mt20=24M
export mt30=36M
export mt40=48M
export mt50=60M
export mt60=72M
export mt70=84M
export mt80=96M
export mt90=108M
export mt100=120M
export hg10=108M
export hg20=96M
export hg30=84M
export hg40=72M
export hg50=60M
export hg60=48M
export hg70=36M
export hg80=24M
export hg90=12M

percent_total

echo "-- Produce 140M total pooled read sets --\n"

export TOTAL=140M

export mt10=14M
export mt20=28M
export mt30=42M
export mt40=56M
export mt50=70M
export mt60=84M
export mt70=98M
export mt80=112M
export mt90=126M
export mt100=140M
export hg10=126M
export hg20=112M
export hg30=98M
export hg40=84M
export hg50=70M
export hg60=56M
export hg70=42M
export hg80=28M
export hg90=14M

percent_total

echo "-- Produce 160M total pooled read sets --\n"

export TOTAL=160M

export mt10=16M
export mt20=32M
export mt30=48M
export mt40=64M
export mt50=80M
export mt60=96M
export mt70=112M
export mt80=128M
export mt90=144M
export mt100=160M
export hg10=144M
export hg20=128M
export hg30=112M
export hg40=96M
export hg50=80M
export hg60=64M
export hg70=48M
export hg80=32M
export hg90=16M

percent_total

echo "-- Produce 180M total pooled read sets --\n"

export TOTAL=180M

export mt10=18M
export mt20=36M
export mt30=54M
export mt40=72M
export mt50=90M
export mt60=108M
export mt70=126M
export mt80=144M
export mt90=162M
export mt100=180M
export hg10=162M
export hg20=144M
export hg30=126M
export hg40=108M
export hg50=90M
export hg60=72M
export hg70=54M
export hg80=36M
export hg90=18M

percent_total

echo "-- Produce 200M total pooled read sets --\n"

export TOTAL=200M

export mt0=0M
export mt10=20M
export mt20=40M
export mt30=60M
export mt40=80M
export mt50=100M
export mt60=120M
export mt70=140M
export mt80=160M
export mt90=180M
export mt100=200M
export hg0=200M
export hg10=180M
export hg20=160M
export hg30=140M
export hg40=120M
export hg50=100M
export hg60=80M
export hg70=60M
export hg80=40M
export hg90=20M
export hg100=0M

percent_total