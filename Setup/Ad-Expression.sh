#!/bin/bash -l

# Author: Diana Trujillo ditrujillo@gmail.com

module load bowtie/2.3.0
module load sratoolkit/2.8.2
module load tophat/2.0.13
module load samtools/1.3
module load cufflinks/2.0.0
module load fastqc/0.11.5
module load bioperl/1.6.901


#Download Arachis duranensis genome
mkdir -p ~/LSE_pipeline/Setup/Ad/fastqc
cd ~/LSE_pipeline/Setup/Ad
wget https://peanutbase.org/files/genomes/Arachis_duranensis/assembly/Aradu_v1.0_20140908.fa.gz
gzip -d Aradu_v1.0_20140908.fa.gz
#Convert to 60 characters per line
perl ../../scripts/60Convert2Fa.pl Aradu_v1.0_20140908.fa
mv Aradu_v1.0_20140908.fa.60 Ad.fa; rm Aradu_v1.0_20140908.fa


#Build Bowtie Index
bowtie2-build Ad.fa Ad


#Download data
for i in SRR2135589 SRR2135593 SRR2135597 SRR2135557 SRR2135560 SRR2135539 SRR2135540 SRR2135601; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR213/$i/$i.sra
 fastq-dump $i.sra --split-3; rm $i.sra
 fastqc $i\_1.fastq; mv *fastqc.* fastqc/.
 fastqc $i\_2.fastq; mv *fastqc.* fastqc/.
done
for i in SRR2140727; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR214/$i/$i.sra
 fastq-dump $i.sra --split-3; rm $i.sra
 fastqc $i\_1.fastq; mv *fastqc.* fastqc/.
 fastqc $i\_2.fastq; mv *fastqc.* fastqc/.
done


#Tophat
#Fastqc reports look fine, so proceed without trimming
###Nodule
cat SRR2135589_1.fastq SRR2135593_1.fastq SRR2135597_1.fastq > AdNods_1.fastq
cat SRR2135589_2.fastq SRR2135593_2.fastq SRR2135597_2.fastq > AdNods_2.fastq
tophat -o Nodules -i 20 -I 2000 --b2-very-sensitive -p 16 Ad AdNods_1.fastq AdNods_2.fastq
###Root
cat SRR2135557_1.fastq SRR2135560_1.fastq > AdRoots_1.fastq
cat SRR2135557_2.fastq SRR2135560_2.fastq > AdRoots_2.fastq
tophat -o Root -i 20 -I 2000 --b2-very-sensitive -p 16 Ad AdRoots_1.fastq AdRoots_2.fastq
###Leaf
cat SRR2135539_1.fastq SRR2135540_1.fastq > AdLeaf_1.fastq
cat SRR2135539_2.fastq SRR2135540_2.fastq > AdLeaf_2.fastq
tophat -o Leaf -i 20 -I 2000 --b2-very-sensitive -p 16 Ad AdLeaf_1.fastq AdLeaf_2.fastq
###Flower
cat SRR2135601_1.fastq SRR2140727_1.fastq > AdFlower_1.fastq
cat SRR2135601_2.fastq SRR2140727_2.fastq > AdFlower_2.fastq
tophat -o Flower -i 20 -I 2000 --b2-very-sensitive -p 16 Ad AdFlower_1.fastq AdFlower_2.fastq

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



