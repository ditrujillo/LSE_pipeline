#!/bin/bash -l

module load fastqc/0.11.5
module load sratoolkit/2.8.2
module load bowtie/2.3.0
module load tophat/2.0.13
module load samtools/1.3
module load cufflinks/2.0.0
module load ncbi_blast+/2.2.28
module load mcl/10.201 
module load bioperl/1.6.901
module load kent/2.4.7
module load clustalw/2.1

#########################################################################
################ Stage2: Expand Pool of candidate genes #################
#########################################################################

################################################
################ Stage2: Expand2 ###############
################################################

################## MCL Cluster #################

#Import data for MCL
mkdir -p ~/LSE_pipeline/Stage2/Expand2/MCL ~/LSE_pipeline/Stage2/Expand2/Blast
cd ~/LSE_pipeline/Stage2/Expand2/MCL
cp ~/LSE_pipeline/Stage2/Expand1/Spada/FilterSpadaHits/SpadaHitsNodSizeMature.fa goodProteins.fasta
cp ~/LSE_pipeline/Stage2/Expand1/Spada/FilterSpadaHits/SpadaHitsNodSize.fa SpadaHitsNodSize.fa

#Blast All-vs.-All set of candidate genes
cd ../Blast
makeblastdb -dbtype prot -in ../MCL/goodProteins.fasta -out goodProteins.index
blastp -db goodProteins.index -query ../MCL/goodProteins.fasta -out MCLblast.tbl -evalue 1e-1 -outfmt 7 -num_threads 16

#Convert Blast Results to MCL format
cd ../MCL
cp -s  ../Blast/MCLblast.tbl .
perl ../../../scripts/mcxdeblast --m9 MCLblast.tbl
mcxassemble -b MCLblast.tbl -r max --map

#Run mcl then obtain labels for the cluster file
#i=1.5 will create bigger clusters
#i=1.8 clusters more tightly
mcl MCLblast.tbl.sym -I 1.8 -o MCLblast.mcl
mcxdump -icl MCLblast.mcl -tabr MCLblast.tbl.tab -o MCLblast.mcl.dump
mkdir ../Groups
orthomclMclToGroups annotated 1 < MCLblast.mcl.dump > ../Groups/groups.txt

################### Find LSEs ##################

#Count ocurrences within each group, then get 'Total Expansion' estimate for each group
cd ../Groups
perl -n -E '$Gm = () = /Gm/g; $Pv = () = /Pv/g; $Mt = () = /Mt/g; $Tp = () = /Tp/g; $Ad = () = /Ad/g; 
            $Tot=$Gm+$Pv+$Mt+$Tp+$Ad; 
            say "Group$.\t$Gm\t$Pv\t$Mt\t$Tp\t$Ad\t$Tot"; ' groups.txt > counts.txt
awk -v OFS="\t" 'BEGIN{print "ID\tGm\tPv\tMt\tTp\tAd\tTot\tRelExp\tIDs"}
   {min=$2; max=$2; for(i=2;i<=(NF-1);++i){min=(min<$i)?min:$i; max=(max>$i)?max:$i} 
   if(min==0){min=1}
   print $0,max/min}' counts.txt > counts-TotExp.txt

#Obtain groups with >=3 total expansion (67 groups)
awk '$8>=3{print $1}' counts-TotExp.txt > ExpansionList.txt

#Get putative function of grouped proteins
#Get IDs & protein sequences of grouped proteins
orthomclExtractProteinIdsFromGroupsFile groups.txt > groups.ProtIDs.txt
perl ../../../scripts/ExtractFastaRecords.pl ../MCL/SpadaHitsNodSize.fa groups.ProtIDs.txt > groups.Prot.fa

#Annotate unknowns with Medicago genes: ~/Data/Mt4.0/Data/genome/proteins.fa
blastp -query groups.Prot.fa -db ~/Data/Mt4.0/Data/genome/MtProtFunction.db -evalue 1e-1 -out Blast.tbl -outfmt 6 -num_threads 16 -max_target_seqs 10
sort -k1,1 -k12,12gr -k11,11g -k3,3gr Blast.tbl | sort -u -k1,1 --merge > bestHits.tbl
awk 'NR==FNR {h[$1] = $2; next} {print $1,h[$1]}' bestHits.tbl groups.ProtIDs.txt > function.txt
perl ../../../scripts/Transpose.pl counts.txt function.txt > countsPlusIDs.txt
awk -v OFS="\t" '{for(i=7; i<=NF; ++i){Func[$i]+=1};{printf $1 "\t"; for (k in Func){printf "("Func[k]")"k", "};print ""; delete Func}}' countsPlusIDs.txt > Summary-IDs.txt
awk -v FS="\t" 'NR==FNR {h[$1] = $2; next} {print $0 "\t"h[$1]}' Summary-IDs.txt counts-TotExp.txt > Counts-TotExp-IDs.txt

################# SPADA search #################

#Align Group Fasta files with ClustalW
mkdir -p ../Spada/Seqs ../Spada/Aln ; cd ../Spada/Seqs
perl ../../../../scripts/SplitFastaByNum.pl ../../Groups/groups.Prot.fa ../../Groups/counts.txt 1000
cd ../Aln
for i in $(cat ../../Groups/ExpansionList.txt); do 
  clustalw2 -infile=../Seqs/$i.fa -type=protein -outfile=$i.aln
done
rm -rf ../Seqs

#Create SPADA HMM profiles
mkdir -p ../SPADAprofile/aln/; cd ../SPADAprofile/aln/
cp -s ../../Aln/Group*.aln .
#cd ~/software/spada/
perl build_profile.pl --cfg ~/LSE_pipeline/ExampleFiles/Cluster.txt --aln ~/LSE_pipeline/Stage2/Expand2/Spada/SPADAprofile/aln --hmm ~/LSE_pipeline/Stage2/Expand2/Spada/SPADAprofile/

#Run Spada
for i in Tp Mt Pv Gm Ad; do
 perl spada.pl --cfg /home/youngn/truji033/LSE_pipeline/ExampleFiles/Cluster.txt \
   --dir /home/youngn/truji033/LSE_pipeline/Stage2/Expand2/Spada/$i \
   --hmm /home/youngn/truji033/LSE_pipeline/Stage2/Expand2/Spada/SPADAprofile/ \
   --fas /home/youngn/truji033/LSE_pipeline/Setup/$i/$i.fa --org $i --evalue 10
done

rm -rf ~/LSE_pipeline/Stage2/Expand2/Spada/*/11_motif_mining ~/LSE_pipeline/Stage2/Expand2/Spada/*/21_model_prediction

########### Expression of Spada Hits ###########

mkdir -p ~/LSE_pipeline/Stage2/Expand2/Spada/Expression
cd ~/LSE_pipeline/Stage2/Expand2/Spada/Expression

#Copy files to run Cufflinks on
#Find tissue expression for the SPADA hits
for j in Tp Mt Pv Gm Ad; do
  mkdir -p $j; cd  $j
  cp -s ../../$j/31_model_evaluation/61_final.g* .
  for i in Flower Root Leaf Nodules; do 
    cufflinks -u --min-intron-length 20 --max-intron-length 2000 -G 61_final.gff -p 16 -o $i \
       --max-bundle-frags 50000000 ~/LSE_pipeline/Setup/$j/$i/accepted_hits.sorted.bam; 
  done; cd ..
done

############ Nod-specific SPADA hits ###########

#Unite Spada IDs to fpkm values for each tissue
for i in Gm Mt Ad Tp Pv; do for j in Flower Leaf Root Nodules; do
  awk -F "\t" 'NR==FNR {FPKM[$1]=$10; next} {if (FNR==1) {next}; gsub("*", "", $19); h[$17]+=1 ; \
       print $17"-"h[$17],$19,$2,FPKM[$2]}' $i/$j/genes.fpkm_tracking $i/61_final.gtb > $i/$j-summary.txt
done; done

for i in Gm Mt Ad Tp Pv; do 
  cd $i
  let "num = 1"; let "next = 2"
  for j in Flower Leaf Root Nodules; do
    awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule"}{print $1}' $j-summary.txt > summary.txt.ALL-1
    awk -v OFS="\t" 'NR==FNR {h[$1]=$4; next} {print $0,h[$1]}' $j-summary.txt summary.txt.ALL-$num > summary.txt.ALL-$next
    let "num = num+1"; let "next = next+1"
  done
  mv summary.txt.ALL-$num ../$i-summary.txt.ALL; rm summary.txt.ALL-*
  cd ..
done

#Get Fasta file with all Spada Hits (4802 peptides)
mkdir -p ../FilterSpadaHits
cd ../FilterSpadaHits
for i in Ad Gm Mt Pv Tp; do
  awk -v OFS="\t" -v Sp=$i '{print ">"Sp"|"$1"\n"$2}' \
                   ../Expression/$i/Flower-summary.txt >> SpadaHits.fa
done

#Get Nod-Specific genes (1434 peptides) with size between 35-250 aa (1400 peptides)
cd ../Expression/
awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule\tAvgNon\tRelExp\tSpecificity"}' > summary.txt.ALL
for i in Gm Mt Ad Tp Pv; do 
    awk -F "\t" -v Sp=$i '{if ($1~"Group") {Spec=0; avgnon=($2+$3+$4)/3; relexp=$5/(avgnon+0.1)
    if ($5>10 && relexp > 4) {Spec=1}
    print Sp"|"$0, avgnon, relexp, Spec}}' $i-summary.txt.ALL >> summary.txt.ALL
done
awk '$8==1{print $1}' summary.txt.ALL > NodSpecList.txt

cd ../FilterSpadaHits
perl ../../../../scripts/ExtractFastaRecords.pl SpadaHits.fa ../Expression/NodSpecList.txt > SpadaHitsNod.fa
faFilter -minSize=35 -maxSize=250 SpadaHitsNod.fa SpadaHitsNodSize.fa		   

#Get mature protein after signal peptide cleavage (1042 peptides)
perl ../../../../scripts/splitJobSigP.pl -p signalp4 -n 1000 SpadaHitsNodSize.fa
cat mature* > SpadaHitsNodSizeMature.fa
perl ../../../../scripts/Collapse100.pl SpadaHitsNodSizeMature.fa
awk -v FS=" " '{print $1}' SpadaHitsNodSizeMature.fa.uniq > SpadaHitsNodSizeMature.fa

