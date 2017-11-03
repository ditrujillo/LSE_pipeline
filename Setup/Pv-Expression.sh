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
mkdir -p ~/LSE_pipeline/Setup/Pv/fastqc
cd ~/LSE_pipeline/Setup/Pv
#Connect to JGI Phytozome 
#https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome#
#Download to ~/LSE_pipeline/Setup/Pv: Pvulgaris_442_v2.0.fa.gz
gzip -d Pvulgaris_442_v2.0.fa.gz
mv Pvulgaris_442_v2.0.fa Pv.fa


#Build Bowtie Index
bowtie2-build Pv.fa Pv


#Download data
for i in SRR1569488 SRR1569478 SRR1569464 SRR1569274; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR156/$i/$i.sra
 fastq-dump $i.sra; rm $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
done


#Tophat
#Fastqc reports look fine, so proceed without trimming 
#Segment length set to 18 and only 1 mismatch allowed due to small read length
#Nodules 21dpi: SRR1569488
tophat -o Nodules -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Pv SRR1569488.fastq
#Root 21dpi: SRR1569478
tophat -o Root -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Pv SRR1569478.fastq
#Flower: SRR1569464
tophat -o Flower -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Pv SRR1569464.fastq
#Leaf 21dpi: SRR1569274
tophat -o Leaf -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Pv SRR1569274.fastq


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


