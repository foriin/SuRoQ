#!/bin/env bash

PS4='+ $(date "+%s.%N")\011 '
set -x

# Function to display usage
usage() {
    echo "Usage: $0 <reads_file> <genome_fasta> <transposon_annotation_fasta> [NCORES] [OUTPUT_DIR]"
    exit 1
}

# Check if correct number of arguments provided
if [ "$#" -lt 3 ]; then
    usage
fi

# Assign arguments to variables
READS_FILE="$1"
GENOME_FASTA="$2"
TRANSPOSON_FASTA="$3"
NCORES="${4:-1}"
SUROQ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
OUT_DIR="${5:-./}"
READS_BASE="$( basename "${READS_FILE}" | cut -f 1 -d '.' )"

# Create output folder if it doesn't exist
mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/idx"
mkdir -p "${OUT_DIR}/map"
mkdir -p "${OUT_DIR}/tables"
mkdir -p "${OUT_DIR}/plots"

# Step 2: Convert reads to .insert format
# This step also removes small RNAs that have homopolimer stretches
echo "Converting reads to .insert format..."
${SUROQ_DIR}/bin/fastx_to_insert $READS_FILE ${OUT_DIR}/${READS_BASE}_reads.insert

# Step 3: Prepare Bowtie index for genome and map reads
if [[ ! -f "${OUT_DIR}/idx/genome_index.1.ebwt" ]]; then
	echo "Preparing Bowtie index for genome..."
	bowtie-build --threads ${NCORES} $GENOME_FASTA ${OUT_DIR}/idx/genome_index 1>/dev/null
else
	echo "Using existing Bowtie index for genome."
fi

echo "Mapping reads to genome..."
bowtie -v 0 -r -a -p ${NCORES} -x ${OUT_DIR}/idx/genome_index -r ${OUT_DIR}/${READS_BASE}_reads.insert --al ${OUT_DIR}/map/${READS_BASE}_reads_genome.insert -S 2>${OUT_DIR}/map/${READS_BASE}_genome.map.log | samtools sort -@ ${NCORES} - -o ${OUT_DIR}/map/${READS_BASE}_reads_genome.bam

echo "Calculating percentage of reads mapped to genome..."
MAPPED_READS_GENOME=$(grep -m 1 "reads with at least one alignment:" ${OUT_DIR}/map/${READS_BASE}_genome.map.log | awk '{print $9}' | tr -d '()' )
echo -e "GENOME\t${MAPPED_READS_GENOME}" > ${OUT_DIR}/tables/${READS_BASE}_map_stats.log

# Prepare size distribution of mapped reads
echo "Preparing size distribution of mapped reads..."
# bedtools bamtobed -i ${OUT_DIR}/map/${READS_BASE}_reads_genome.bam > ${OUT_DIR}/map/${READS_BASE}_reads_genome.bed
# ${SUROQ_DIR}/bin/bed_to_bed2 ${OUT_DIR}/${READS_BASE}_reads.insert ${OUT_DIR}/map/${READS_BASE}_reads_genome.bed > ${OUT_DIR}/map/${READS_BASE}_reads_genome.bed2
# sort -k7,7 ${OUT_DIR}/map/${READS_BASE}_reads_genome.bed2 | cut -f 7,4 | uniq > ${OUT_DIR}/map/${READS_BASE}_reads_genome.insert
awk '{print length($1)}' ${OUT_DIR}/map/${READS_BASE}_reads_genome.insert | sort | uniq -c > ${OUT_DIR}/tables/${READS_BASE}_genome_uq_sd.txt
awk '{len=length($1); freq=$2; for(i=0;i<freq;i++) print len}' ${OUT_DIR}/map/${READS_BASE}_reads_genome.insert | sort -n | awk '{count[$1]++} END {for (len in count) print count[len], len}' | sort -k2,2n > ${OUT_DIR}/tables/${READS_BASE}_genome_all_sd.txt

# Step 5: Prepare Bowtie index for TEs and map reads
if [[ ! -f "${OUT_DIR}/idx/te_index.1.ebwt" ]]; then
	echo "Preparing Bowtie index for transposable elements (TEs)..."
	bowtie-build --threads ${NCORES} $TRANSPOSON_FASTA ${OUT_DIR}/idx/te_index 1>/dev/null
else
	echo "Using existing Bowtie index for TEs."
fi

echo "Mapping reads to TEs..."
bowtie -v 0 -r -a -p ${NCORES} -x ${OUT_DIR}/idx/te_index -r ${OUT_DIR}/${READS_BASE}_reads.insert -S 2>${OUT_DIR}/map/${READS_BASE}_te.map.log | samtools sort - -o ${OUT_DIR}/map/${READS_BASE}_reads_te.bam
echo "Calculating percentage of reads mapped to TEs..."
MAPPED_READS_TE=$(grep -m 1 "reads with at least one alignment:" ${OUT_DIR}/map/${READS_BASE}_te.map.log | awk '{print $9}' | tr -d '()' )
echo -e "TE\t${MAPPED_READS_TE}" >> ${OUT_DIR}/tables/${READS_BASE}_map_stats.log

# Convert SAM to BED and BED2
bedtools bamtobed -i ${OUT_DIR}/map/${READS_BASE}_reads_te.bam > ${OUT_DIR}/map/${READS_BASE}_reads_te.bed
${SUROQ_DIR}/bin/bed_to_bed2 ${OUT_DIR}/${READS_BASE}_reads.insert ${OUT_DIR}/map/${READS_BASE}_reads_te.bed > ${OUT_DIR}/map/${READS_BASE}_reads_te.bed2
# Get ping-pong signature data
echo "Getting ping-pong signature data..."
${SUROQ_DIR}/bin/ping_pong -a ${OUT_DIR}/map/${READS_BASE}_reads_te.bed2 -b ${OUT_DIR}/map/${READS_BASE}_reads_te.bed2 -p ${NCORES} > ${OUT_DIR}/tables/${READS_BASE}_te_mapped.pp
# Get Size Distro and PFM for + strand
echo "Getting size distributions and PFMs for TE sense-mapped reads... (only unique reads)"
# samtools view -F 20 ${OUT_DIR}/map/${READS_BASE}_reads_te.bam | awk '{print length($10)}' | sort | uniq -c > ${OUT_DIR}/tables/${READS_BASE}_te_sense_sd.txt
awk '$6 == "+" {print $7 "\t" $4}' ${OUT_DIR}/map/${READS_BASE}_reads_te.bed2 | sort | uniq > ${OUT_DIR}/map/${READS_BASE}_reads_te_sense.insert
awk '{print length($1)}' ${OUT_DIR}/map/${READS_BASE}_reads_te_sense.insert | sort | uniq -c > ${OUT_DIR}/tables/${READS_BASE}_te_sense_sd.txt
${SUROQ_DIR}/bin/insert_to_pfm ${OUT_DIR}/map/${READS_BASE}_reads_te_sense.insert > ${OUT_DIR}/tables/${READS_BASE}_te_map_sense.pfm

# Get Size Distro and PFM for - strand
echo "Getting size distributions and PFMs for TE antisense-mapped reads... (only unique reads)"
# samtools view -f 16 ${OUT_DIR}/map/${READS_BASE}_reads_te.bam | awk '{print length($10)}' | sort | uniq -c > ${OUT_DIR}/tables/${READS_BASE}_te_asense_sd.txt
awk '$6 == "-" {print $7 "\t" $4}' ${OUT_DIR}/map/${READS_BASE}_reads_te.bed2 | sort | uniq > ${OUT_DIR}/map/${READS_BASE}_reads_te_asense.insert
awk '{print length($1)}' ${OUT_DIR}/map/${READS_BASE}_reads_te_asense.insert | sort | uniq -c > ${OUT_DIR}/tables/${READS_BASE}_te_asense_sd.txt
${SUROQ_DIR}/bin/insert_to_pfm ${OUT_DIR}/map/${READS_BASE}_reads_te_asense.insert > ${OUT_DIR}/tables/${READS_BASE}_te_map_asense.pfm

# Prepare a plot
echo "Drawing plots..."
Rscript ${SUROQ_DIR}/bin/draw_plots.R ${OUT_DIR} ${READS_BASE}

# echo "${OUT_DIR} ${READS_BASE}"

echo "Kudos for a job... done."

 


