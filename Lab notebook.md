# Ecological Genomics Lab notebook    

## Author: Aayudh Das     

### We are going to discuss various aspects of Transcriptomics

### Date started: (2017-01-18)   
### Date end:   (2017-05-02)

### Table of contents    
* [Page 1: 2017-01-18](#id-section1). Ecological genomics, first class
* [Page 2: 2017-01-23](#id-section2). The next-generation sequencing 
* [Page 3:2017-01-25](#id-section3). Thinking critically about genomic discovery of ecological adaptation
* [Page 4:2017-01-30](#id-section4). Group presentations of project ideas
* [Page 5:2017-02-01](#id-section5). Sequencing strategies applied to biological questions: WGS, RNA-Seq; RAD/GBS, Amplicons
* [Page 6:2017-02-06](#id-section6). Transcriptomics 1
* [Page 7:2017-02-08](#id-section7). Transcriptomics 2
* [Page 8:2017-02-13](#id-section8). Transcriptomics 3
* [Page 9:2017-02-15](#id-section9). Transcriptomics 4
* [Page 10:2017-02-20](#id-section10).Transcriptomics 5
* [Page 11:2017-02-22](#id-section11). Transcriptomics 6
* [Page 12:2017-03-01](#id-section12). Catch up day and WGCNA package
* [Page 13:2017-03-06](#id-section13).Population genomics 1
* [Page 14:2017-03-08](#id-section14).Population genomics 2 and ASSIGNMENT 2 codes
* [Page 14:2017-03-20](#id-section15).Population genomics 3

------
<div id='id-section1'/>
### Page 1: 2017-01-18. Ecological genomics, first class

### **Steve and Melissa's intro**    
* Introduction to transcriptomics. The goals of learning in this course. Lots of RNA seq analysis.
* Ecological genomics institute, KSU: emphasis on adaptation to environment   
* Gordon Research Conference: Integrating different levels of biological organization on any system. approach and tool focused! Field going towards new data and new analytic techniques  
* Intro to eco genomics, oxford press; Using technology to address ecological issues such as nutrient cycling, population structure, life history vairation , trophic interaction, stress responess, and adpatation to environmental change   

------
<div id='id-section2'/>
### Page 2: 2017-01-23. The next-generation sequencing

### **Discussion format**
* Outline-20min time, Engaging activity, Use board effectively, Example from literature.

### Info update: Next generation sequencing 

#### Advances in sequence Tech
* Human Genome project 2001-2003, Sanger, 15 years, 1genome, $3B
* AiSeqXTen-2014 Illumina, 1day, 45 whole genomes in a day, $100 each.

#### Range of Application
* Whole genome sequencing (WGS), RNA seq, ChIP seq, targeted/ capture seq (design a probe) 
* Where is genetic variation Phenotypes?
* Number of samples 
* Comparative studies already modeled or not (1st thing to consider)
* Demographic activity, adaptive genetic variation, gene expression, length of the reads, 
* Number of reads (short 50bp, long 100, 150, 300 bp in MiSeq 10000-60000 SMRT; reads single vs paired end), distribution

#### General Library Prep workflow
* Extraction of DNA or RNA to cDNA
* Fragment sample
* Ligate adaptors
* Individual barcode
* Add seq adaptors
* PCR

#### Sequencing by synthesis (SBS)
* Bridge amp, cluster generation, labeled dNTP (ATCG)

#### Other technologies

#### Learning activity

### Discussion paper: Ellegren 2014

```
[Ellegren et al. 2014] (http://www.sciencedirect.com.ezproxy.uvm.edu/science/article/pii/S0169534713002310)
```
Go to the paper [Ellegren 2014] (http://www.sciencedirect.com.ezproxy.uvm.edu/science/article/pii/S0169534713002310). 

* Genome sequencing and population genomics in non-model organisms
* High-throughput sequencing technologies are revolutionizing the life sciences.

------
<div id='id-section3'/>
### Page 3: 2017-01-25. Thinking critically about genomic discovery of ecological adaptation

### Info update: Quantitative trait nucleotides (QTNs)

#### What are QTN?
Identifying quantitative trait nucleotide; example trait-flowering time, flower color, venom potency, tolerance, defense compounds, toxin tol, drought tol, altitude tol.  

#### Quantitative genetic theory of adaptive traits 
* Va
* h2

#### Methods
* Linkage mapping
* GWAS
* Selection scan

### Discussion paper: Rockman 2012; Lee et al.2014

```
[Rockman 2012] (https://www-ncbi-nlm-nih-gov.ezproxy.uvm.edu/pmc/articles/PMC3386609/)
```
Go to the paper [Rockman 2012] (https://www-ncbi-nlm-nih-gov.ezproxy.uvm.edu/pmc/articles/PMC3386609/)

```
[Lee et al.2014] (https://academic-oup-com.ezproxy.uvm.edu/aobpla/article/doi/10.1093/aobpla/plu004/156429/Identifying-the-genes-underlying-quantitative)
```
Go to the paper [Lee et al.2014] (https://academic-oup-com.ezproxy.uvm.edu/aobpla/article/doi/10.1093/aobpla/plu004/156429/Identifying-the-genes-underlying-quantitative)

------
<div id='id-section4'/>
### Page 4: 2017-01-30. Group presentations of project ideas

### Discussion on Sea Star RNA seq data set
* one page proposal due in 2nd Feb

#### Our grop idea:
* Immune related gene expression
  * Reverse pathology
  * Looking specific classes of genes
  * A prior test of resistance genes
      * Compare individual that stayed healthy vs those got sick.
      * Looking at the sick to healthy transition.

##### We are looking the gene expression, up arms paper.
* H vs S (day 9, 12, 15)
* H vs H (day 9, 12, 15)
* S vs S (day 9, 12, 15) 

```
[Project proposal] (https://docs.google.com/document/d/13YSQoPu9I8blGWvmKFSttaG6CcuCs3YCZ7mMc9_hxxc/edit?ts=58928c12#heading=h.59c9l1pt0wvm)
```
Propsal link [Project proposal] (https://docs.google.com/document/d/13YSQoPu9I8blGWvmKFSttaG6CcuCs3YCZ7mMc9_hxxc/edit?ts=58928c12#heading=h.59c9l1pt0wvm)

------
<div id='id-section5'/>
### Page 5: 2017-02-01. Sequencing strategies applied to biological questions: WGS, RNA-Seq; RAD/GBS, Amplicons

### INFO UPDATE

It was me who was giving all the updates

### Terminal outputs
* How to login to the terminal?
  * aayudhdas$ ssh aadas@pbio381.uvm.edu
  * Then give your UVM password to login

* Command: top 
  * Shows you whoever is working in the server
  * to quit: [aadas@pbio381 ~]$ q

* To see your present working directory
  * [aadas@pbio381 ~]$ pwd
* To see what's the current contents of any folders and files in our current location
  * [aadas@pbio381 ~]$ ll

* Let’s make a new folder (aka, directory) using the mkdir command. Let’s name this folder “mydata”
```  
[aadas@pbio381 ~]$ mkdir mydata
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ mkdir scripts
```
* We can change our current location within the directory structure using the cd command. Let’s use cd to move inside the mydata/ directory and ll to list its contents:
``` 
[aadas@pbio381 ~]$ cd mydata
[aadas@pbio381 mydata]$ ll
total 0
[aadas@pbio381 mydata]$ cd /data/
[aadas@pbio381 data]$ ll
```
* We’ve placed the text file containing all the metadata information on the seastar sampling under a shared space on the server. The path to this shared space is:
  * /data/ Try using cd to navigate over to this location. Then ll to show its contents. You should see something like this:
``` 
[aadas@pbio381 mydata]$ cd /data/
[aadas@pbio381 data]$ ll
```

* Now, cd into the folder called “project_data” and ll. Do you see this?
``` 
[aadas@pbio381 data]$ cd project_data/
[aadas@pbio381 project_data]$ ll
```
* The file called “ssw_samples.txt” is the one with the seastar metadata. We don’t want to open and make changes to this file in the shared space, because we don’t want to have our edits affect the rest of the group. So, let’s first make a copy of this file over to our home directory and put it inside the “mydata” folder. Use the cp command, followed by the filename, and the path to your destination (remember the ~ signals your home directory, and each subdirectory is then separated by a /):
``` 
[aadas@pbio381 project_data]$ cp ssw_samples.txt ~/mydata
```
* cd back to your ~/mydata/ directory and look inside. You should see your file… 
``` 
[aadas@pbio381 project_data]$ cd ~/mydata/
[aadas@pbio381 mydata]$ ll
```
* Let’s take a peek at this file with the head command, which prints the first 10 lines to screen.
``` 
[aadas@pbio381 mydata]$ head ssw_samples.txt 
Customize your lines-
[aadas@pbio381 mydata]$ head -n 15 ssw_samples.txt 
```
* Let’s see the last 10 lines
  * [aadas@pbio381 mydata]$ tail ssw_samples.txt 

* What if we want to extract just the rows of data that correspond to Healthy (HH) individuals? We can use the search tool grep to search for a target query. Any line matching our search string will be printed to screen.
``` 
[aadas@pbio381 mydata]$ grep 'HH'ssw_samples.txt
[aadas@pbio381 mydata]$ grep 'SS' ssw_samples.txt
```
* What if instead of printing it to screen, we want to save the output of our search to a new file? This is easy, just use the “>” symbol to redirect the results of any command to an output file with your choice of name.
``` 
[aadas@pbio381 mydata]$ grep 'SS' ssw_samples.txt >ssw_SSonly.txt
[aadas@pbio381 mydata]$ ll
[aadas@pbio381 mydata]$ mkdir sample_by_disease/
[aadas@pbio381 mydata]$ grep 'HH' ssw_samples.txt >ssw_HHonly.txt
[aadas@pbio381 mydata]$ ll
```

------
<div id='id-section6'/>
### Page 6: 2017-02-06. Transcriptomics 1

### INFO UPDATE


### Terminal codes

```
#how to login
ssh aadas@pbio381.uvm.edu
And then give your netID pwd

# how to go to the directory
[aadas@pbio381 ~]$ cd /data/project_data/fastq

# do ls to view all the file
[aadas@pbio381 fastq]$ ll

#dataset I will be working
-rw-r--r--. 1 mlloyd   users  756575578 Feb  2 12:21 10_5-08_H_0_R1.fq.gz
-rw-r--r--. 1 mlloyd   users  820537378 Feb  2 12:21 10_5-08_H_0_R2.fq.gz

# To view the file
[aadas@pbio381 fastq]$ zcat 10_5-08_H_0_R1.fq.gz | head

#navigate to 
[aadas@pbio381 fastq]$ cd /data/scripts/
[aadas@pbio381 scripts]$ ll

# Just copied 
[aadas@pbio381 scripts]$ cp trim_example.sh ~/scripts/

[aadas@pbio381 scripts]$ cd ~
[aadas@pbio381 ~]$ ll
total 0
drwxr-xr-x. 3 aadas users 93 Feb  1 11:29 mydata
drwxr-xr-x. 2 aadas users 36 Feb  6 11:07 scripts
[aadas@pbio381 ~]$ cd scripts/
[aadas@pbio381 scripts]$ ll
total 4
-rwxr--r--. 1 aadas users 711 Feb  6 11:07 trim_example.sh

#to see the file
[aadas@pbio381 scripts]$ head trim_example.sh 

#vim and file name to edit it
vim trim_example.sh 

#edit the script
    /data/project_data/fastq/10_5-08_H_0_R1.fq.gz \
                 /data/project_data/fastq/10_5-08_H_0_R2.fq.gz \
                /data/project_data/fastq/cleanreads/"samp_R1_clean_paired.fq" \
                /data/project_data/fastq/cleanreads/"samp_R1_clean_unpaired.fq" \
                /data/project_data/fastq/cleanreads/"samp_R2_clean_paired.fq" \
                /data/project_data/fastq/cleanreads/"samp_R2_clean_unpaired.fq" \

#to save
i for insert
Esc then :w to save
q: for quit
:wq! Save quit at the same time

#how to run the file
[aadas@pbio381 scripts]$ bash trim_example.sh 
# open a new window and login and do top to see whether it’s running or not.

#running stat
Input Read Pairs: 13271979 Both Surviving: 11213144 (84.49%) Forward Only Surviving: 1510836 (11.38%) Reverse Only Surviving: 262066 (1.97%) Dropped: 285933 (2.15%)

# to rename the file
mv samp_R1_clean_unpaired.fq 10_5-08_H_0_R1_clean_unpaired.fq 

# to create the html file
[aadas@pbio381 cleanreads]$ fastqc 10_5-08_H_0_R
10_5-08_H_0_R1_clean_paired.fq    10_5-08_H_0_R2_cleaned_paired.fq  
10_5-08_H_0_R1_clean_unpaired.fq  10_5-08_H_0_R2_clean_unpaired.fq  
[aadas@pbio381 cleanreads]$ fastqc 10_5-08_H_0_R
10_5-08_H_0_R1_clean_paired.fq    10_5-08_H_0_R2_cleaned_paired.fq  
10_5-08_H_0_R1_clean_unpaired.fq  10_5-08_H_0_R2_clean_unpaired.fq  
[aadas@pbio381 cleanreads]$ fastqc 10_5-08_H_0_R*

#to copy the html file to my desktop
aayudhdas$ scp aadas@pbio381.uvm.edu:/data/project_data/fastq/cleanreads/10_5-08_H_0_R1_clean_paired_fastqc.html ~/Desktop/

# only the paired files should be considered.

# copied one file from the class directory to my folder
cp command 
```
------
<div id='id-section7'/>
### Page 7: 2017-02-08. Transcriptomics 2

### INFO UPDATE

#### Introduction
* Why wild systems: 
     * Non-model and non-traditional model organisms.
     * Silent genes responding to multiple stimuli
     * Novel transcript without homologs and closely related model organism

#### Brief overview of transcriptomics
*	Microarray 
     * Easy for ecological analysis
       *RNA seq: 
     * genome wide ecological transcripts, more involved analysis.

#### Main questions
* Variation is there in gene expression and how’s it structured?
     * Evolutionary process
     * Gene expression is heritable (natural selection)
     * Epigenetics 
     * Qst-Fst comparison
     * qQTL expression and quantitative trait loci mapping 

* How does the gene expression affect phenotype?

* How to environmental stimuli affect gene expression?
     * Abiotic stress
     * Env heterogeneity 
     * Host parasite interaction
     * Selective biotic and abiotic interaction 
     * Looking at the molecular level
          * Molecular level
          * Genotype
          * Phenotypic plasticity
     * Some limitation while working in the environment 
          * You need to flash freeze your sample
          * You are only getting a snapshot

* How does gene expression affect phenotype?
     * Alternative phenotype
     * Moving from transgenic, RNA and CRISPER/CAS

#### Problems 
* Bias in signal
* Polyploidy
* RNA pooling
* Statistical analysis
* Unannotated genes

------
<div id='id-section8'/>
### Page 8: 2017-02-13. Transcriptomics 3

### INFO UPDATE

#### Background 
* Enables Dt examination (inter poplulation individual)
  * Disease resistance 
  * Mating behavious 
  * Adaptive behaviour

* Molecular Mechanism 
  * Phenotypical/ behaviroal plasticity (migration pattern)

#### Limitation
* Refernce genome quality gene annotation availablility expense per sample library preparation

#### Issues
* Under utilization of biological replicates
  * Requiring indep. library preparations
  * Doesn't involved pooled samples
  * 23/158 studies 
  * Derive broad bio colclusions

* Priotonize sequence depth over replication (problem)

* Nice dynamic range of RNA seq data-- Noisy
  * Poisson counting error
  * Biological variance
  * Technical variance

#### Rules of thumb
* Use of more biological replicates instead of depth
* Seq. depth > 10 reads / transcript 
   * ~10-20 Million mapped reads/sample
* 3 Biological replicates per condition
* Conduct a pilot experiment- 
   * What's the best/ powerful experiment I can afford? 
   * What's the smallest fold change I conduct?

### Discussion paper: Johnston et al. 2016

```
[Johnston et al. 2016 (http://onlinelibrary.wiley.com/doi/10.1111/mec.13879/abstract)
```
Go to the paper [Johnston et al. 2016] (http://onlinelibrary.wiley.com/doi/10.1111/mec.13879/abstract). 

------
<div id='id-section9'/>
### Page 9: 2017-02-15. Transcriptomics 4

### INFO UPDATE: SNPs and population genomics
* SNP data-expressed sequences

#### Process
```
Tissue--> Sequence-->Clean/trim-->assembly-->SNP detection-->validation
```
##### Tissue
* breadth of tissue, developmental stages, exon sapping

##### Pool & sequencing libraries
* ~30-100 million paired ending read

##### Clean trim  
* Process raw seq. data
* Important for SNP detection
* Digital normalization
  * Remove high coverage reads and associated errors
  * loss of quantitative information
* Assemble cleaned paired long reads
* Prune-reduce DNA contaminaiton, non-coding RNA, gene fragments

##### Assembly
* Evaluation-COGS

##### SNP detection
* Software-constant patterns of sequence validation
  * Sequence error: Eliminate SNPs of low frequency
  * Artifacts caused by InDels
  * Quality score

##### SNP validation 
* Designing primers
* Sequencing by Mass spectrometry 

#### Application
* Difference in population structure
* How natural selection acting on particular loci

##### Methods of application 
* Outliers: For a given locus what's the level of differentiation comapred to diff. across genome.
* Non-outliers: Tests high FST loci for other factors associated with selection.
  * Fitness
  * Functional enrichment

### Discussion paper: Zhao et al. 2016

```
[Zhao et al. 2016 (https://academic.oup.com/mbe/article/33/3/707/2579453/Global-Transcriptional-Profiling-of-Diapause-and)
```
Go to the paper [Zhao et al. 2016] (https://academic.oup.com/mbe/article/33/3/707/2579453/Global-Transcriptional-Profiling-of-Diapause-and). 

### Terminal Codes (15th Feb):

How to go to the script folder

```
ip0af52760:~ aayudhdas$ ssh aadas@pbio381.uvm.edu
[aadas@pbio381 ~]$ cd scripts/

[aadas@pbio381 scripts]$ ll
total 8830312
-rw-r--r--. 1 aadas users 8802028313 Feb 13 11:33 10_5-20_S_2_bwaaln.sam
-rw-r--r--. 1 aadas users  119846624 Feb 13 11:20 10_5-20_S_2_R1.fq.gz_left_clean_paired.fq.sai
-rw-r--r--. 1 aadas users  120349760 Feb 13 11:30 10_5-20_S_2_R2.fq.gz_right_clean_paired.fq.sai
-rwxr-xr-x. 1 aadas users        902 Feb 13 11:09 bwaaln.sh
-rwxr--r--. 1 aadas users        765 Feb  6 11:26 trim_example.sh
```

Let’s check out our .sam files!

```
[aadas@pbio381 scripts]$ tail -n 100 10_5-20_S_2_bwaaln.sam > tail.sam
[aadas@pbio381 scripts]$ vim tail.sam
:set nowrap
```

Let’s see how many of our reads map uniquely

```
[aadas@pbio381 scripts]$ grep -c XT:A:U 10_5-20_S_2_bwaaln.sam
3036521
[aadas@pbio381 scripts]$ grep -c X0:i:1 10_5-20_S_2_bwaaln.sam
3056482
```

to go back to my directory

```
[aadas@pbio381 scripts]$ cd ~/
```

Extract read counts from the .sam file from each sample

```
[aadas@pbio381 scripts]$ cd /data/scripts/
[aadas@pbio381 scripts]$ cp countxpression_PE.py ~/scripts/
```

To open it in python

```
[aadas@pbio381 scripts]$ sed -i 's/::/|_/g' 10_5-20_S_2_bwaaln.sam
[aadas@pbio381 scripts]$ python countxpression_PE.py 20 35 countstatssummary.txt 10_5-20_S_2_bwaaln.sam

#things to remeber 
"sed" find and replace command
"-i" find and replace right in the command line
"s" search
"/::/" search for two colons
"/ \ _ /" replace with an underscore
g = option
```

see my file

```
[aadas@pbio381 scripts]$ ll
total 8830952
-rw-r--r--. 1 aadas users     611483 Feb 15 11:33 10_5-20_S_2_bwaaln_counts.txt
-rw-r--r--. 1 aadas users 8802028313 Feb 15 11:27 10_5-20_S_2_bwaaln.sam
-rw-r--r--. 1 aadas users  119846624 Feb 13 11:20 10_5-20_S_2_R1.fq.gz_left_clean_paired.fq.sai
-rw-r--r--. 1 aadas users  120349760 Feb 13 11:30 10_5-20_S_2_R2.fq.gz_right_clean_paired.fq.sai
-rwxr-xr-x. 1 aadas users        902 Feb 13 11:09 bwaaln.sh
-rw-r--r--. 1 aadas users        224 Feb 15 11:33 countstatssummary.txt
-rwxr--r--. 1 aadas users       6487 Feb 15 11:13 countxpression_PE.py
-rw-r--r--. 1 aadas users      27269 Feb 15 10:37 tail.sam
-rwxr--r--. 1 aadas users        765 Feb  6 11:26 trim_example.sh
```

------

<div id='id-section10'/>

### Page 10: 2017-02-22. Transcriptomics 5

### INFO UPDATE: none 

### Terminal code

First go to your current directory 

``` 
[aadas@pbio381 ~]$ pwd
/users/a/a/aadas
[aadas@pbio381 ~]$ ll
```

Let's go to my data

Copy  all the files from one directory to my location. "." signifies to the location that you are presently in.

```
[aadas@pbio381 mydata]$ scp /data/project_data/DGE/* .
[aadas@pbio381 mydata]$ ll
```
Now move to your MAC

Open a new terminal, now that's your hard drive. I want to copy everything to my own folder in MAC

```
ip0af52d66:~ aayudhdas$ pwd
/Users/aayudhdas
ip0af52d66:~ aayudhdas$ cd ~/Dropbox/Aayudh_UVM/ecological\ genomics/
```

Now lets make a folder in my dropbox

```	
ip0af52d66:ecological genomics aayudhdas$ mkdir RNA_seq
ip0af52d66:ecological genomics aayudhdas$ cd RNA_seq
ip0af52d66:RNA_seq aayudhdas$ ls		
```

Now copy everything to the new folder

```
ip0af52d66:RNA_seq aayudhdas$ scp aadas@pbio381.uvm.edu:/users/a/a/aadas/mydata/* .
aadas@pbio381.uvm.edu's password: 
cols_data_trim.txt                            100% 1982     1.0MB/s   00:00    
countsdata_trim.txt                           100% 5420KB   9.9MB/s   00:00    
countstatsummary.txt                          100% 8122     2.4MB/s   00:00    
DESeq2_exploreSSW_trim.R                      100% 9018     2.7MB/s   00:00    
explore_expression_data.R                     100%  439   204.4KB/s   00:00    
scp: /users/a/a/aadas/mydata/sample_by_disease: not a regular file
samples_by_disease                            100%  462    24.3KB/s   00:00    
ssw_samples.txt                               100% 1255   223.0KB/s   00:00    
ip0af52d66:RNA_seq aayudhdas$ 
```

## R Script

```
source("http://bioconductor.org/workflows.R")
workflowInstall("rnaseqGene")

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

setwd("~/Dropbox/Aayudh_UVM/ecological genomics/RNA_seq")
library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#################### Build dataset, model, and run analyses

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ day + location + health)
# In this typical model, the "sex effect" represents the overall effect controlling for differences due to population and devstage. page 26-27 manual DESeq2.pdf
# The last term in the model is what is tested.  In this case sex.
#  This is not the same as an interaction.


dim(dds)
[1] 26550    65

dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
[1] 13334    65  # at > 100; little more than an average of 10 reads per sample for the 93 samples

colSums(counts(dds))
hist(colSums(counts(dds)), breaks = 80, xlim=c(0,max(colSums(counts(dds)))))


colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S"))

dds <- DESeq(dds)  # this step takes a loooong time ~4 minutes with the trimmed data set
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3308 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

save(dds, file="dds.trim.Robject")

res <- results(dds)
res <- res[order(res$padj),]
head(res)
# log2 fold change (MAP): health S vs H 
# Wald test p-value: health S vs H 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
# TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i1_g.23138_m.23138 1950.0719       2.488783 0.4311875
# TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110  902.2693       2.475891 0.4599085
# TRINITY_DN43359_c0_g1_TRINITY_DN43359_c0_g1_i1_g.14658_m.14658  889.9707       1.163219 0.2482335
# TRINITY_DN47215_c1_g4_TRINITY_DN47215_c1_g4_i3_g.25054_m.25054  774.1126       1.723917 0.3650258
# TRINITY_DN47215_c0_g1_TRINITY_DN47215_c0_g1_i5_g.25051_m.25051  911.7634       1.586693 0.3431307
# TRINITY_DN45416_c4_g2_TRINITY_DN45416_c4_g2_i3_g.19333_m.19333 1629.8753       1.775765 0.3873817
# stat       pvalue         padj
# <numeric>    <numeric>    <numeric>
# TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i1_g.23138_m.23138  5.771927 7.837024e-09 8.691260e-06
# TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110  5.383443 7.307426e-08 4.051968e-05
# TRINITY_DN43359_c0_g1_TRINITY_DN43359_c0_g1_i1_g.14658_m.14658  4.685987 2.786136e-06 7.724563e-04
# TRINITY_DN47215_c1_g4_TRINITY_DN47215_c1_g4_i3_g.25054_m.25054  4.722727 2.327027e-06 7.724563e-04
# TRINITY_DN47215_c0_g1_TRINITY_DN47215_c0_g1_i5_g.25051_m.25051  4.624166 3.761091e-06 8.342101e-04
# TRINITY_DN45416_c4_g2_TRINITY_DN45416_c4_g2_i3_g.19333_m.19333  4.584018 4.561241e-06 8.430694e-04

summary(res)
# out of 13334 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 50, 0.37% 
# LFC < 0 (down)   : 8, 0.06% 
# outliers [1]     : 539, 4% 
# low counts [2]   : 11686, 88% 
# (mean count < 80)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

plotMA(res, main="DESeq2", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i1_g.23138_m.23138", intgroup=(c("status","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = date)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + scale_y_log10(breaks=c(25,100,1000)) + ylim(0,2500)
p

# 2.2 Data quality assessment by sample clustering and visualization 

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))

# rld <- rlog(dds, blind=FALSE) # this takes too long with such a large data set!
# plotPCA(rld, intgroup=c("status","date"))
```

------

<div id='id-section11'/>

### Page 11: 2017-02-27. Transcriptomics 6

### INFO UPDATE:  

None 

### Discussion paper: 

1. Coalescnt 
2. reticulation
3. Purigying/background selection
4. Gene trees vs species 
5. Introgression / recombination 
6. Incomplete linkage

### Terminal Code

First open your own UVM server and do this

```
[aadas@pbio381 ~]$ scp /data/project_data/DGE/* .
cp: omitting directory ‘/data/project_data/DGE/round1’
[aadas@pbio381 ~]$ ll
```

Now open a new window

```
ip0af5257c:~ aayudhdas$ pwd
/Users/aayudhdas
ip0af5257c:~ aayudhdas$ cd ~/Dropbox/Aayudh_UVM/ecological\ genomics/
ip0af5257c:ecological genomics aayudhdas$ mkdir RNA_seq1
ip0af5257c:ecological genomics aayudhdas$ cd RNA_seq1
ip0af5257c:RNA_seq1 aayudhdas$ ls
ip0af5257c:RNA_seq1 aayudhdas$ scp aadas@pbio381.uvm.edu:/data/project_data/DGE/* .
aadas@pbio381.uvm.edu's password: 
allcountsdata.txt                                                100% 3883KB   9.6MB/s   00:00    
cols_data_trim.txt                                               100% 2342   893.1KB/s   00:00    
countsdata_trim2.txt                                             100% 3394KB  12.1MB/s   00:00    
countstatsummary.txt                                             100% 8059     3.0MB/s   00:00    
DESeq2_SSW_round2.R                                              100%   13KB   4.0MB/s   00:00    
scp: /data/project_data/DGE/round1: not a regular file
ip0af5257c:RNA_seq1 aayudhdas$ 
```
------
<div id='id-section12'/>

### Page 12: 2017-03-01. Catch up day and 'R' package update

### INFO UPDATE: WGCNA (weighted gene correlation network analysis) package  

#### Outline

##### Overview of WGCNA

Describe correlation (co-expression) patterns among genes in micro-array samples.

Network construction —> Module identification —> Relationship of modules to external info—> Relationship between/within modules—> Finding key drivers in modules of interest 

##### Network construction

Node=Gene

Edge=Strength of co-rrleation in expression

Signed network-positive correlation in expression

Unsigned network- Absolute value of correlation in expression (+ve and -ve correlation)

##### Module detection

Unspecified clustering

1. Incorporation of excel info into network
2. Topological properties
3. Other features
4. Limitation

### Data anlysis

LTR= log likelihood in R

It's going to take liklihood of 4 models and then it will take the likelihood of reduced model.

LR=chi square, df=#pair reduced

response variable: Gene expression

#### For practice, let's work with fewer genes so that the models can run faster...

dds <- dds[sample(nrow(dds), 1200), ] comment out that 

-----

<div id='id-section13'/>

### Page 13: 2017-03-06. Population genomics

### INFO UPDATE: Population genetics

#### Population genomics: 

##### SNPs

##### sampling unit is individuals 

#### Processes: 

1. Population structure 
2. Diversity with populations
3. Selections 
   1. Positive
   2. Negative; "purifying"

#### Pipeline 

Raw reads—>clean—> Assemble —> Mapped reads—>Count reads—>DGE

​											    —> Call SNPs + Genotypes (population genomics)

##### Challenges of SNP coding

1.  Sequencing error (Illumina 1:100)

##### Filters

1. Minor allele frequency 

   #### Challenges of calling genotypes (AA, AT, TT) 

   1. Multinomial distribution
   2. Predicts A=T=0.5
   3. 100% hetero (SNPs where its filtered out)

#### pi heterozygocity

For sequence of i and j, pi = Xi XjPiiy

Pisyn = 4Neu

Pinonsyn/Pisyn = measure of how strong the selection is

### Terminal Code

```
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps/
[aadas@pbio381 reads2snps]$ ll
[aadas@pbio381 reads2snps]$ ll *vcf
[aadas@pbio381 reads2snps]$ vim head_SSW_bamlist.txt.vcf 
:set nowrap
[aadas@pbio381 reads2snps]$ vcftools --vcf SSW_bamlist.txt.vcf 
[aadas@pbio381 reads2snps]$ grep "unres" SSW_bamlist.txt.vcf | wc
5631864 185851488 1028494934
[aadas@pbio381 reads2snps]$ grep "para" SSW_bamlist.txt.vcf | wc
   4354  143652  795592
[aadas@pbio381 reads2snps]$ vcftools --vcf SSW_bamlist.txt.vcf --min-alleles 2 --max-alleles 2
[aadas@pbio381 reads2snps]$ vcftools --vcf SSW_bamlist.txt.vcf --maf 0.02
[aadas@pbio381 reads2snps]$ vcftools --vcf SSW_bamlist.txt.vcf --max-missing 0.8
```

wc= word count

Now do it all this together

```
[aadas@pbio381 reads2snps]$ vcftools --vcf SSW_bamlist.txt.vcf --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/biallelic.MAF0.02.Miss0.8
```

now go to your home directory

```
[aadas@pbio381 data]$ cd ~/
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ vim biallelic.MAF0.2.Miss0.8.recode.vcf
:set nowrap
```

To see the mannual 

```
[aadas@pbio381 ~]$ man vcftools
```

to hardy

```
[aadas@pbio381 ~]$ vcftools --vcf biallelic.MAF0.2.Miss0.8.recode.vcf --hardy
[aadas@pbio381 ~]$ head out.hwe
```

Use R from this

```
[aadas@pbio381 ~]$ R
> getwd()
[1] "/data/users/a/a/aadas"
> hardy <-read.table("out.hwe", header=T) 
> str(hardy)
```

p value

```
> hardy[which(hardy$P_HET_EXCESS<0.001),]
> hardy[which(hardy$P_HET_DEFICIT<0.001),]
>hardy[which(hardy$P_HET_DEFICIT<0.001),2:7]
```

to quit R

```
> quit()
Save workspace image? [y/n/c]: y
```

```
[aadas@pbio381 ~]$ vcftools --vcf biallelic.MAF0.2.Miss0.8.recode.vcf --geno-r2
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ vim out.geno.ld 
[aadas@pbio381 ~]$ R

> LD <- read.table("out.geno.ld", header=T)
> str(LD)
'data.frame':	4650 obs. of  5 variables:
 $ CHR   : Factor w/ 192 levels "TRINITY_DN35598_c0_g1_TRINITY_DN35598_c0_g1_i1_g.5802_m.5802",..: 105 105 105 105 105 105 105 105 105 105 ...
 $ POS1  : int  2127 2127 2127 2127 2127 2127 2127 2127 2127 2127 ...
 $ POS2  : int  2217 2235 2244 2276 2277 2535 2805 2970 2994 3327 ...
 $ N_INDV: int  19 19 19 19 19 20 19 20 20 20 ...
 $ R.2   : num  0.00654 0.00309 0.00309 0.00309 0.00309 ...
> LD$dist<- abs(LD$POS1-LD$POS2)
> str(LD)
'data.frame':	4650 obs. of  6 variables:
 $ CHR   : Factor w/ 192 levels "TRINITY_DN35598_c0_g1_TRINITY_DN35598_c0_g1_i1_g.5802_m.5802",..: 105 105 105 105 105 105 105 105 105 105 ...
 $ POS1  : int  2127 2127 2127 2127 2127 2127 2127 2127 2127 2127 ...
 $ POS2  : int  2217 2235 2244 2276 2277 2535 2805 2970 2994 3327 ...
 $ N_INDV: int  19 19 19 19 19 20 19 20 20 20 ...
 $ R.2   : num  0.00654 0.00309 0.00309 0.00309 0.00309 ...
 $ dist  : int  90 108 117 149 150 408 678 843 867 1200 ...
> pdf("LD.plot.pdf")
> plot(LD$dist,LD$R.2)
> dev.off()

```
------
<div id='id-section14'/>
### Page 14: 2017-03-08. Population genomics

### INFO UPDATE: Population genomics

#### Rate of evolution due to relationship

##### Ne=Effective popultion size

4 methods:

1.  From species life history
2.  From varioance in allele frequency between generation
3.  Genetic polymorphism
4.  Correlated trait ~ body size 

Ne vary across:

1.  Species
2.  Genome
    1. Genetic hitchhiking~selective slope
    2. Background selection
    3. Fewer sex chromosomes if outcross goes high.

##### Mutation

1. Across the whole gene and chromosome 
   1. Duplication
   2. Inversion
   3. Delection
   4. Translocation
2. Mutation at the base level i.e. point mutation (main source is substitution)
   1. Transition (Purine to Purine and Purine to pyriidine)
   2. Transversion (purone to purine)

Now base level mutation are two types

1.  Synonymous (silent)
2.  Non-synonymous (Change in the DNA) / Natural selection operates
    1. Purifying selection
    2. Positive selection  

5 classes

1. Neutral mutation: Effect of fitness<1 
2. slightly deletorious 
3. Deletorius 
4. Advantages

##### Relationship between Ne and Substitution rate (NeRR)

1. Selective sweeps
2. Clonal interference 

##### Take home:

1. The study of NeRR help us to understand the process that drives and limits evolution


2. Drift and selection are the most important forces determining NeRR.
3. We need to work on better way to estimate Ne.
4. With time in DNA sed will let us estimate better Ne, substitute rate, mutation rates.
5. Most mutations are deletorious. 

#### Terminal Code

```
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps/
[aadas@pbio381 reads2snps]$ ll
[aadas@pbio381 reads2snps]$ vcftools --gzvcf SSW_byind.txt.vcf.gz 
```

After filtering, kept 22 out of 22 Individuals

After filtering, kept 7485987 out of a possible 7485987 Sites

Run Time = 20.00 seconds

```
[aadas@pbio381 reads2snps]$ vcftools --gzvcf SSW_byind.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.8
```

After filtering, kept 5565 out of a possible 7485987 Sites

```
[aadas@pbio381 reads2snps]$ cd ~/
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ gzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --hardy
```

After filtering, kept 5565 out of a possible 5565 Sites

```
[aadas@pbio381 ~]$ head out.hwe
```

in R now

```
> hwe<-read.table("out.hwe", header=T)
> str(hwe)
> summary(hwe)
> which(hwe$P_HET_DEFICIT<0.01)
```

[1] 1001 1021 1023 1300 1302 1320 1407 1409

```
> hwe[which(hwe$P_HET_DEFICIT<0.01),]
```

read the file

```
[aadas@pbio381 ~]$ vi ssw_healthloc.txt
then press o
[aadas@pbio381 reads2snps]$ cat ssw_healthloc.txt 
Individual	Trajectory	Location	SNPs
10	HH	INT	N
24	HH	INT	Y
27	HH	INT	Y
08	HS	INT	Y
09	HS	INT	Y
15	HS	INT	Y
19	HS	INT	Y
20	HS	INT	Y
03	SS	INT	Y
07	SS	INT	Y
14	SS	INT	Y
22	SS	INT	Y
23	SS	INT	Y
26	SS	INT	Y
28	SS	INT	Y
29	SS	INT	Y
31	HH	SUB	Y
32	HH	SUB	Y
33	HH	SUB	Y
34	HH	SUB	N
35	HH	SUB	Y
36	SS	SUB	Y
37	MM	SUB	Y
38	MM	SUB	Y

```

```
[aadas@pbio381 reads2snps]$ grep "HH" ssw_healthloc.txt > ~/H_OneSampPerInd.txt
[aadas@pbio381 reads2snps]$ grep "SS" ssw_healthloc.txt > ~/S_OneSampPerInd.txt
[aadas@pbio381 reads2snps]$ grep "HS" ssw_healthloc.txt >> ~/S_OneSampPerInd.txt
```

go to directory

```
[aadas@pbio381 ~]$ cut -f 1 H_OneSampPerInd.txt >H_OneSampPerInd2.txt 
[aadas@pbio381 ~]$ cat H_OneSampPerInd2.txt
10
24
27
31
32
33
34
35
[aadas@pbio381 ~]$ cut -f 1 S_OneSampPerInd.txt >S_OneSampPerInd2.txt 
[aadas@pbio381 ~]$ cat S_OneSampPerInd2.txt
03
07
14
22
23
26
28
29
36
08
09
15
19
20
```

To see the allele frequency

```
[aadas@pbio381 ~]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz  --freq2 --keep H_OneSampPerInd2.txt --out H_AlleleFreqs
```

After filtering, kept 5565 out of a possible 5565 Sites

```
[aadas@pbio381 ~]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz  --freq2 --keep S_OneSampPerInd2.txt --out S_AlleleFreqs
```

After filtering, kept 5565 out of a possible 5565 Sites

```
[aadas@pbio381 ~]$ vim H_AlleleFreqs.frq
:set nowrap
edit by i
then esc
then :wq to save 
```

### Homework codes

For 'int'

```
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#Subsetting colData
i=which(colData[,5]=="int")
colData1=colData[i,]
#Subsetting countsdataData, below 29 score is "int"
countData1=(countData [,1:48])
#Dimension of the matrices 
dim(countData)
dim(colData)
dim(colData1)
dim(countData1)
#################### MODEL NUMBER 1: TEST EFFECT OF HEALTH CONTROLLING FOR 'int'

dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ + health)
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
dds <- dds[sample(nrow(dds))]
dim(dds)
colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) 
dds <- DESeq(dds) 
res <- results(dds)
res <- res[order(res$padj),]
head(res)
summary(res)

################################################################
# Data quality assessment by sample clustering and visualization 

plotMA(res, main="DESeq2", ylim=c(-3,3))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,1000)
p


## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
p

############## PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))
```

for 'sub'

```
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#Subsetting colData
i=which(colData[,5]=="sub")
colData2=colData[i,]
#Subsetting countsdataData, 49-77 column is "sub"
countData2=(countData [,49:77])
#Dimension of the matrices 
dim(colData)
dim(countData)
dim(colData2)
dim(countData2)
#################### MODEL NUMBER 1: TEST EFFECT OF HEALTH CONTROLLING FOR LOCATION

dds <- DESeqDataSetFromMatrix(countData = countData2, colData = colData2, design = ~ + health)
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
dds <- dds[sample(nrow(dds))]
dim(dds)
colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) 
dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)
summary(res)


################################################################
# Data quality assessment by sample clustering and visualization 

plotMA(res, main="'sub' location", ylim=c(-6,5))

## Check out one of the genes to see if it's behaving as expected
## change the ylim based on the scale
d <- plotCounts(dds, gene="TRINITY_DN29480_c0_g1_TRINITY_DN29480_c0_g1_i1_g.3557_m.3557", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,100)
p


## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN29480_c0_g1_TRINITY_DN29480_c0_g1_i1_g.3557_m.3557", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
p

############## PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))
```

both

```
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#Dimension of the matrices 
dim(colData)
dim(countData)

#################### MODEL NUMBER 1: TEST EFFECT OF HEALTH CONTROLLING FOR LOCATION

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
dds <- dds[sample(nrow(dds))]
dim(dds)
colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) 
dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)
summary(res)


################################################################
# Data quality assessment by sample clustering and visualization 

plotMA(res, main="both location", ylim=c(-3,3))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,5000)
p


## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
p

############## PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))
```


=======
-----
<div id='id-section15'/>
### Page 15: 2017-03-20. Population genetics 3

### INFO UPDATE: Population genetics 3

### Introduction

#### Model datastructure

admixture: Maximmum likelihood approach

### Methods for GA (global ancestry) estimation

#### Non parametric approach

1. **Clustering**: pairise matrix; need some program
2. PCA, multi-scale dimensioning (MSD) 

### Method of LA (local ancestry) estimation

1. Hidden-Markov model
2. **Program**: LAMP; Alternatively RFMix-uses discrimination of functionl analysis. 
3. It's very accurate for 2 ways mixture.


### Weakness 

we know the no. of population and allele frequency. But we don't actually know. To fix this we need to simulate.  

### Future perspective 

Genomics will produce dense SNPs.

