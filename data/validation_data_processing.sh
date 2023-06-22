#!/bin/bash

# Check RNA-Seq
if [ $# -eq 0 ]; then
    echo "Please provide as the first argument the path to the RNA-Seq FASTQ file"
    exit 2
fi
FASTQ_FILE=$1
HUMAN_RNA_REF=GRcH38_human.fna
THREADS=`grep -c ^processor /proc/cpuinfo`

# Remove FASTQ information from the RNA-Seq FASTQ file
echo 'Convert FASTQ to FASTA, keeping only the sequence'
cat ${FASTQ_FILE} | sed -n --expression='1~4s/^@/>/p;2~4p' > RNASeq.fasta

# Prepare database from the human transcriptome GRcH38 (NCBI)
if [ ! -f "$HUMAN_RNA_REF" ]; then
    echo "Reference human transcriptome was not found, download it"
    wget -O- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz" | gunzip > $HUMAN_RNA_REF
fi

if [ ! -f "HUMANRNA.ndb" ]; then
    makeblastdb -in $HUMAN_RNA_REF -input_type fasta -dbtype nucl -title human -parse_seqids -out HUMANRNA
fi

# Search against the genome
echo 'Start BLASTN search, this needs some time... 8 hours with 16 threads for a 4.5 GB FASTQ file'
blastn -query RNASeq.fasta -num_threads $THREADS -db HUMANRNA -evalue 0.01 -perc_identity 70 -outfmt "6 qseqid sseqid qstart qend sstart send nident pident evalue bitscore score mismatch gapopen gaps qcovs qcovhsp slen qlen length" | gzip -c >  HUMAN_RNA_RESULT.txt.gz

# Sorting results # pident evalue
echo 'Start sorting results'
zcat HUMAN_RNA_RESULT.txt.gz | sort --version-sort -k1,1 -k8,8 -k9,9r | sort --version-sort -u -k1,1 | gzip -c > HUMAN_RNA_RESULT_SORTED.txt.gz

# Remove redundant matches, ensure unique sequences
echo 'Remove redundant matches'
zcat HUMAN_RNA_RESULT_SORTED.txt.gz | sort --version-sort -u -k2,2 | gzip -c > HUMAN_RNA_RESULT_SORTED_REDUNDANCY.txt.gz

# Reconstruction the transcript sequences from the matches with the coding and non-coding label from the NCB reference transcriptome data
echo 'Execute MOSGA ClassifySeqFromBlast module'
if [ ! -f "mosga" ]; then
    git clone --depth 1 https://gitlab.com/mosga/mosga.git
fi
python3 mosga/accumulator/mosga.py misc ClassifySeqFromBlast HUMAN_RNA_RESULT_SORTED_REDUNDANCY.txt.gz,GRcH38_human.fna,RNASeq.fasta,./ -vvvv

# The output will be NexusDNAGSK_SC_coding_redundancy_removed.fasta and NexusDNAGSK_SC_noncoding_redundancy_removed.fasta, respectively coding.fast and noncoding.fasta
