# OConnell_2018

##

**Prerequisites:**  
Cutadapt 1.19  
TopHat 2.1.1  
Bowtie 2.2.7  
featureCounts 1.6


Human genomic sequences and annotation files (GRCh38.p12) were downloaded from the [NCBI repository](ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/).  

| files             | MD5 check sum (unzipped)         | Description                                               |
| ----------------- |:--------------------------------:| ----------------------------------------------------------|
| GRCh38.p12.fa     | a4cac7d7ac4dd31ac68b384b10cf444d | RNA in fasta format, coding + noncoding                   |
| GRCh38.p12.fna    | 860290186a4ee3e95cd48dc528a45363 | Genome sequence, chromosomes and extrachromosomal contigs |
| GRCh38.p12.gbk    | 3c35b07e638485984479d50dd5cfebca | RNA in gene bank format, coding + noncoding               |
| GRCh38.p12.gff    | 56394751c00a5bdfb74152a7ed146855 | Genome annotation                                         |
   
  
**Preparation of genome annotation for gene expression analysis.** Extrachromosomal contigs and annotations were omitted. 'Gnomon' (Predicted) records from gff file were also omitted and only 'RefSeq' (manually curated) left. Perl and R scripts are included in the GitHub repository.   
```bash
Discard_extrachromosomal_contigs.pl GRCh38.p12.fna >GRCh38.p12.custom.fna
Discard_extrachromosomal_annotation.pl GRCh38.p12.gff >GRCh38.p12.custom.gff
Discard_gnomon_annotation.pl >GRCh38.p12.Refseq.gff
```
Non-coding RNA genes were removed, only coding genes with their mRNA, transcript, exon, and CDS children features were left.
```sh
Discard_noncoding_annotation.R
```



**Indexing human genome**
```bash
bowtie2-build ./Human_indices/GRCh38.p12.custom.fna ./Human_indices/NCBI_genome
```
**Indexing mouse transcriptome for TopHat**  
```bash
tophat -G GRCh38.p12.Refseq.coding.gff --transcriptome-index ./tophat-2.1.1/Human_indices/Refseq_coding ./bowtie2-2.2.7/Human_indices/NCBI_genome
```
**Discarding ribosomal, mitochondrial and tRNA reads**
```bash
```
