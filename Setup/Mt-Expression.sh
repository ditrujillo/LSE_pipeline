#!/bin/bash -l

# Author: Diana Trujillo ditrujillo@gmail.com

module load bowtie/2.3.0
module load sratoolkit/2.8.2
module load tophat/2.0.13
module load samtools/1.3
module load cufflinks/2.0.0
module load fastqc/0.11.5
module load bioperl/1.6.901


#Download Medicago truncatula 4.0 genome
mkdir -p ~/LSE_pipeline/Setup/Mt/fastqc
cd ~/LSE_pipeline/Setup/Mt
#Connect to JGI Phytozome 
#http://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Mtruncatula#
#Download to ~/LSE_pipeline/Setup/Mt: Mtruncatula_285_Mt4.0.fa.gz 
gzip -d Mtruncatula_285_Mt4.0.fa.gz 
mv Mtruncatula_285_Mt4.0.fa Mt.fa


#Build Bowtie Index
bowtie2-build Mt.fa Mt


#Download data
for i in SRR350517 SRR350520 SRR350519 SRR350518; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR350/$i/$i.sra
 fastq-dump $i.sra; rm $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
done


#Tophat
#Fastqc reports look fine, so proceed without trimming 
#Segment length set to 18 and only 1 mismatch allowed due to small read length
#Nodules: SRR350517
tophat -o Nodules -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Mt SRR350517.fastq
#Root: SRR350520
tophat -o Root -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Mt SRR350520.fastq
#Flower: SRR350519
tophat -o Flower -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Mt SRR350519.fastq
#Leaf: SRR350518
tophat -o Leaf -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Mt SRR350518.fastq


#Cufflinks
for i in Nodules Root Flower Leaf; do
 cd $i
 samtools sort accepted_hits.bam -o accepted_hits.sorted.bam
 samtools index accepted_hits.sorted.bam
 rm accepted_hits.bam unmapped.bam
 cufflinks -u --min-intron-length 20 --max-intron-length 2000 -p 16 -o Cufflinks --max-bundle-frags 50000000 \
   accepted_hits.sorted.bam
 cd ..
done


