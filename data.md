# Data


## Bed file

Original bed file downloaded from NCBI ftp and filtered for the following genes:
TTN, BRAF, KRAS, OMA1, and TGDS.

```shell
# Download file
wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz

# Process file
gunzip nexterarapidcapture_expandedexome_targetedregions.bed.gz
grep -P "\tTTN|\tBRAF|\tKRAS|\tOMA1|\tTGDS" nexterarapidcapture_expandedexome_targetedregions.bed > gene_panel_exomes.bed
```

## Reads
Fastq files downloaded from NCBI ftp, aligned to GRCh37, and then filtered for the reads that
fall in the genomic regions of the following genes:
TTN, BRAF, KRAS, OMA1, and TGDS.

```shell
# Download fastq files
wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R1_001_trimmed.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R2_001_trimmed.fastq.gz

# Align
# bwa v0.7.17
# samtools v1.9
# main calls:
bwa mem -M -R @RG\tID:NA12878-N1-DNA1-WES1.0\tSM:NA12878-N1-DNA1-WES1\tPL:ILLUMINA -p -t 16 /path/to/GRCh37/hs37d5/hs37d5.fa /dev/stdin
samtools sort -T /path/to/temporary/file/tmp.N7fk2CV2aG/sort_bam -m 4G -@ 4 -O BAM -o work/bwa.NA12878-N1-DNA1-WES1/out/bwa.NA12878-N1-DNA1-WES1.bam

# Merge genomic exome regions
bedtools merge -i gene_panel_exomes.bed -d 100000 > gene_panel_exomes_merged.bed

# Filter bam files and revert to fastq
samtools view -b -h -L gene_panel_exomes.bed bwa.NA12878-N1-DNA1-WES1.bam > gene_panel.bam
bedtools bamtofastq -i gene_panel.bam -fq giab_gene_panel_R1.fastq -fq2 giab_gene_panel_R2.fastq
```
