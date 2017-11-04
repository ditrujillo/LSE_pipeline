LSE Discovery Pipeline
============

CONTENTS
------------

-   Overview
-   Citation
-   Dependencies
-   Usage:
1. Get LSE pipeline scripts
2. Manually download genomes from Phytozome 
3. Run Setup.sh to map reads to genomes, get nodule transcriptome
4. Run Stage1.sh to obtain initial pool of candidate nodule-specific small secreted peptides
5. Run Stage2-Expand1.sh to identify lineage-specific expansions of nodule-specific small secreted peptides
6. Run Stage2-Expand2.sh to refine the peptides identified for candidate LSE families
7. Run Stage3.sh on the command line interactively, in step with the manual curation steps described in LSE_pipeline.sh
8. Run Heatmap.R to combine phylogenetic and expression data to confirm LSEs


OVERVIEW
------------

The **LSE Discovery Pipeline** was developed to detect Lineage Specific Expansions (LSEs) of nodulation-related small secreted peptides across legume species. 

In the first stage, the pipeline uses RNA-seq data from four tissues, including nodules, that is mapped to the reference genome for each species. Transcripts with nodule-enhanced expression are identified, then reading frames for those tissues are filtered for size and presence of a signal peptide. 

In the second stage, the pool of putative nodule-specific secreted peptides are clustered by sequence homology, to identify clusters that are overrepresented in specific legume lineages. Those with lineage-specific expansions are retained, and used to conduct homology searches against the reference genomes. The second stage is repeated twice. In this way, the initial pool of candidate genes is filtered to retain LSEs, and the number of members within those LSE gene families is expanded iteratively.

In the third stage, the candidate LSE families are manually curated through combined phylogenetic and expression analyses.


CITATION
------------

**Author:** Diana Trujillo ditrujillo@gmail.com

"Cross-species examination of legume signaling peptides reveals that nodule-specific PLAT domain proteins are required for nodulation" Trujillo et al. (In Prep.)


DEPENDENCIES
------------

LSE_pipeline.sh is a shell script that requires loading (if working within the University of Minnesota Supercomputing Institue) or installing the following programs (make available in your $PATH):

- fastqc/0.11.5
- sratoolkit/2.8.2
- bowtie/2.3.0
- tophat/2.0.13
- samtools/1.3
- cufflinks/2.0.0
- ncbi_blast+/2.2.28
- mcl/10.201 
- mcxdeblast (from mcl/12.135; download to ~/LSE_pipeline/scripts)
- bioperl/1.6.901
- kent/2.4.7
- clustalw/2.1
- orthomcl/2.0
- SPADA: https://github.com/orionzhou/SPADA
- signalP/4.0


USAGE
------------

**Note:** A detailed explanation of the LSE discovery pipeline steps and required user input is outlined in [LSE_pipeline.sh](LSE_pipeline.sh)

1. Get LSE pipeline scripts

   **cd**
   
   **git clone git://github.com/ditrujillo/LSE_pipeline.git**


2. Manually download genomes from Phytozome 

   **mkdir -p ~/LSE_pipeline/Setup/Mt**
   
   **mkdir -p ~/LSE_pipeline/Setup/Gm**
   
   **mkdir -p ~/LSE_pipeline/Setup/Tp**
   
   **mkdir -p ~/LSE_pipeline/Setup/Pv**

   **Note:** Connect to JGI Phytozome: https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome#

   Download *M. truncatula* genome to ~/LSE_pipeline/Setup/Mt: Mtruncatula_285_Mt4.0.fa.gz 
   
   Download *G. max* genome to        ~/LSE_pipeline/Setup/Gm: Gmax_275_v2.0.fa.gz
   
   Download *T. pratense* genome to   ~/LSE_pipeline/Setup/Tp: Tpratense_385_v2.fa.gz
   
   Download *P. vulgaris* genome to   ~/LSE_pipeline/Setup/Pv: Pvulgaris_442_v2.0.fa.gz


3. Run Setup.sh to map reads to genomes, and obtain genomic boundaries for nodule transcripts

   **cd ~/LSE_pipeline**
   
   **./Setup.sh**


4. Run Stage1.sh to obtain initial pool of candidate nodule-specific small secreted peptides

   **./Stage1.sh**


5. Run Stage2-Expand1.sh to identify lineage-specific expansions of nodule-specific small secreted peptides

   **Note:** This step relies on the Small Peptide Alignment Detection Application (SPADA)

   Prior to running SPADA, the SPADA configuration file ~/LSE_pipeline/ExampleFiles/Cluster.txt must be set up to run on the user's computer, depending on the location of the program dependencies.

   Stage2-Expand1.sh must be modified with the locations of the user files, in order to run spada.pl
 
   **./Stage2-Expand1.sh**


6. Run Stage2-Expand2.sh to refine the peptides identified for candidate LSE families

   **Note:** Stage2-Expand2.sh must be modified with the locations of the user files, in order to run spada.pl

   **./Stage2-Expand2.sh**


7. Run Stage3.sh on the command line interactively, in step with the manual curation steps described in LSE_pipeline.sh

   **Note:** Large gene families are sometimes split into separate clusters after running Stage2-Expand2.sh
   
   The peptide clusters can be grouped into gene families according to #Stage3# in LSE_pipeline.sh

   **./Stage3.sh**


8. Run Heatmap.R to combine phylogenetic and expression data to confirm LSEs

   **./Heatmap.R**




