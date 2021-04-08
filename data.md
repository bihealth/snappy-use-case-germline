# Data


## Bed file

Original bed file downloaded from NCBI ftp and filtered for the following genes:
TTN, BRAF, KRAS, OMA1, and TGDS.

```shell
# Download and filter bed file
wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed
# Process file
grep -P ":TTN|:BRAF|:KRAS|:OMA1|:TGDS" TruSeq_exome_targeted_regions.hg19.bed > gene_panel_exomes.bed
```

## Reads
The aligned BAM file was downloaded and then filtered for the reads that
fall in the genomic regions of the following genes:
TTN, BRAF, KRAS, OMA1, and TGDS.
The resulting file was converted to fastq.

```shell
# Download bam file
wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam

# Merge genomic exome regions
# 1	58946391	59012446
# 2	179390719	179672150
# 7	140433814	140624564
# 12	25358179	25403854
# 13	95226308	95248511
bedtools merge -i gene_panel_exomes.bed -d 100000 > gene_panel_exomes_merged.bed
sed -i "s/chr//g" gene_panel_exomes_merged.bed

# Filter bam files and revert to fastq
bedtools intersect -abam NIST-hg001-7001-ready.bam -b  gene_panel_exomes_merged.bed > gene_panel.bam
samtools index gene_panel.bam
bedtools bamtofastq -i gene_panel.bam -fq giab_gene_panel_R1.fastq -fq2 giab_gene_panel_R2.fastq
```
