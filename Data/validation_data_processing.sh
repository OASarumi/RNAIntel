#!/bin/bash

# Remove FASTQ information from the RNA-Seq NexusDNAGSK_SC-CHE-01_L7.D705_2.clipped.fastq.gz
zcat NexusDNAGSK_SC-CHE-01_L7.D705_2.clipped.fastq.gz | sed -n --expression='1~4s/^@/>/p;2~4p' > NexusDNAGSK_SC-CHE-01_L7.D705_2.clipped.fasta

# Prepare database from the human transcriptome GRcH38 (NCBI)
makeblastdb -in GRcH38_human.fna -input_type fasta -dbtype nucl -title human -parse_seqids -out HUMANRNA

# Search against the genome
blastn -query NexusDNAGSK_SC-CHE-01_L7.D705_2.clipped.fasta -num_threads 15 -db HUMANRNA -evalue 0.01 -perc_identity 70 -outfmt "6 qseqid sseqid qstart qend sstart send nident pident evalue bitscore score mismatch gapopen gaps qcovs qcovhsp slen qlen length" | gzip -c >  HUMAN_RNA_RESULT.txt.gz

# Sorting results # pident evalue
zcat HUMAN_RNA_RESULT.txt.gz | sort --version-sort -k1,1 -k8,8 -k9,9r | sort --version-sort -u -k1,1 | gzip -c > HUMAN_RNA_RESULT_SORTED.txt.gz

# Remove redundant matches, ensure unique sequences
zcat HUMAN_RNA_RESULT_SORTED.txt.gz | sort --version-sort -u -k2,2 | gzip -c > HUMAN_RNA_RESULT_SORTED_REDUNDANCY.txt.gz

# Reconstruction the transcript sequences from the matches with the coding and non-coding label from the NCB reference transcriptome data
python3 /opt/mosga/accumulator/mosga.py misc ClassifySeqFromBlast HUMAN_RNA_RESULT_SORTED.txt.gz,GRcH38_human.fna,NexusDNAGSK_SC-CHE-01_L7.D705_2.clipped.fasta,./ -vvvv

# The output will be NexusDNAGSK_SC_coding_redundancy_removed.fasta and NexusDNAGSK_SC_noncoding_redundancy_removed.fasta, respectively coding.fast and noncoding.fasta
