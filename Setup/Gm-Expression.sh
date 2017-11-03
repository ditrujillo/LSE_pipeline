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
mkdir -p ~/LSE_pipeline/Setup/Gm/fastqc
cd ~/LSE_pipeline/Setup/Gm
#Connect to JGI Phytozome 
#https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome#
#Download to ~/LSE_pipeline/Setup/Gm: Gmax_275_v2.0.fa.gz
gzip -d Gmax_275_v2.0.fa.gz
mv Gmax_275_v2.0.fa Gm.fa


#Build Bowtie Index
bowtie2-build Gm.fa Gm


#Download data
for i in SRR037385 SRR037387 SRR037382 SRR037384; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR037/$i/$i.sra
 fastq-dump $i.sra; rm $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
done
for i in SRR899512 SRR899511 SRR899431 SRR899430; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR899/$i/$i.sra
 fastq-dump $i.sra; rm $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
done


#Tophat
#Fastqc reports look fine, so proceed without trimming 
#Segment length set to 18 and only 1 mismatch allowed due to small read length
#Nodules: SRR037385 SRR899512
cat SRR037385.fastq SRR899512.fastq > Nodules.fastq
tophat -o Nodules -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Gm Nodules.fastq
#Root: SRR037387 SRR899511
cat SRR037387.fastq SRR899511.fastq > Root.fastq
tophat -o Root -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Gm Root.fastq
#Flower: SRR037382 SRR899431
cat SRR037382.fastq SRR899431.fastq > Flower.fastq
tophat -o Flower -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Gm Flower.fastq
#Leaf: SRR037384 SRR899430
cat SRR037384.fastq SRR899430.fastq > Leaf.fastq
tophat -o Leaf -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 18 --segment-mismatches 1 Gm Leaf.fastq

rm SRR213* SRR214* 

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


