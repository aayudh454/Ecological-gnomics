# Ecological Genomics Lab notebook    

## Author: Aayudh Das     

### We are going to discuss various aspects of Transcriptomics

#### Login ID: ssh aadas@pbio381.uvm.edu

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
* [Page 15:2017-03-20](#id-section15).Population genomics 3
* [Page 16:2017-03-22](#id-section16).Population genomics 4
* [Page 17:2017-03-27](#id-section17).Population genomics 5
* [Page 18:2017-03-29](#id-section18).Population genomics 6
* [Page 19:2017-03-31](#id-section19).Homework3_population genetics
* [Page 20:2017-04-03](#id-section20).Catch up day
* [Page 21:2017-04-05](#id-section21).Enrichment and annotation
* [Page 22:2017-04-10](#id-section22).Metagenomics 1
* [Page 23:2017-04-12](#id-section23).Metagenomics 2
* [Page 24:2017-04-17](#id-section24).Metagenomics 3
* [Page 25:2017-04-17](#id-section25).Metagenomics 4
* [Page 25:2017-04-24](#id-section26).Project Code

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

### Terminal codes

```
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps
[aadas@pbio381 reads2snps]$ ll
[aadas@pbio381 reads2snps]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.8
```

After filtering, kept 5317 out of a possible 7486938 Sites

```
[aadas@pbio381 reads2snps]$ cd ~/
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ gzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps/ssw_healthloc.txt
-bash: cd: /data/project_data/snps/reads2snps/ssw_healthloc.txt: Not a directory
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps/ssw_healthloc.txt
-bash: cd: /data/project_data/snps/reads2snps/ssw_healthloc.txt: Not a directory
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps/
[aadas@pbio381 reads2snps]$ ll
[aadas@pbio381 reads2snps]$ cat ssw_healthloc.txt 
[aadas@pbio381 reads2snps]$ grep "HH" ssw_healthloc.txt | cut -f1 >~/H_SampleIDs.txt
[aadas@pbio381 reads2snps]$ cd ~/
[aadas@pbio381 ~]$ cat H_SampleIDs.txt
```

10

24

27

31

32

33

34

35

```
[aadas@pbio381 ~]$ wc H_SampleIDs.txt
 8  8 24 H_SampleIDs.txt
[aadas@pbio381 ~]$ cd /data/project_data/snps/reads2snps/
[aadas@pbio381 reads2snps]$ ll
[aadas@pbio381 reads2snps]$ cat ssw_healthloc.txt
[aadas@pbio381 reads2snps]$ grep "HS\|SS" ssw_healthloc.txt | cut -f1 >~/S_SampleIDs.txt
[aadas@pbio381 reads2snps]$ cd ~/
[aadas@pbio381 ~]$ ll
[aadas@pbio381 ~]$ wc S_SampleIDs.txt
14 14 42 S_SampleIDs.txt
[aadas@pbio381 ~]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep H_SampleIDs.txt --out H_AlleleFreqs
```

After filtering, kept 5317 out of a possible 5317 Sites

```
[aadas@pbio381 ~]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep S_SampleIDs.txt --out S_AlleleFreqs
```

After filtering, kept 5317 out of a possible 5317 Sites

```
vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --weir-fst-pop H_SampleIDs.txt --weir-fst-pop S_SampleIDs.txt --out HvS_Fst
```

After filtering, kept 5317 out of a possible 5317 Sites

#### Transfer the file

```
#open a new window in your mac
ip0af5268c:users aayudhdas$ cd ~/Dropbox/Aayudh_UVM/ecological\ genomics/
# create new folder
ip0af5268c:ecological genomics aayudhdas$ mkdir pop_gen
#now transfer the file from the location, your directory is=pwd
ip0af5268c:pop_gen aayudhdas$ scp aadas@pbio381.uvm.edu:/users/a/a/aadas/H_AlleleFreqs.frq .
```

R part

```
setwd("/Users/aayudhdas/Dropbox/Aayudh_UVM/ecological genomics/pop_gen")
list.files()

H_freq <- read.table("H_AlleleFreqs.frq", header=T)
S_freq <- read.table("S_AlleleFreqs.frq", header=T)

str(H_freq)
All_freq <- merge(H_freq, S_freq, by=c("CHROM", "POS"))
str(All_freq)
head(All_freq)
All_freq$diff <- (All_freq$H_ALT - All_freq$S_ALT)
hist(All_freq$diff, breaks=50, col="red", main = "Allele frequency difference (H-S)")

fst <- read.table("HvS_Fst.weir.fst", header=T)
All_freq.fst <- merge(All_freq, fst, by=c("CHROM", "POS"))
plot(All_freq.fst$diff, All_freq.fst$WEIR_AND_COCKERHAM_FST, xlab="Allele frequency difference (H-S)", ylab="Fst", main="Healthy vs. Sick SNP divergence")
# Which are the genes that are showing the highest divergence between Healthy and Sick?
All_freq.fst[which(All_freq.fst$WEIR_AND_COCKERHAM_FST>0.2),]

All_freq.fst[which(All_freq.fst$CHROM=="TRINITY_DN47139_c0_g1_TRINITY_DN47139_c0_g1_i1_g.24693_m.24693"),]

points( 0.250000,0.2098250,bg="red", cex=2)
```

## =======

<div id='id-section16'/>

### Page 16: 2017-03-22. Population genetics 4

### INFO UPDATE: Species divergence with gene flow

#### 1. Allopatric vs sympatric

**Allopatric speciation**-species divergence in absence of gene flow; physical barrier present

**Sympatric speciation** is the process through which new species evolve from a single ancestral species while inhabiting the same geographic region; Gene flow present; selected genes appear diverging; related alleles appear homogenous. 

#### 2. Inferring history of divergence 

1.   **Genomic scans:** 

     look island of differentiation (distribution, sumamry statistics that easures differentiation. High FST values region under selection) 

     -gene vs population trees (compare assumed population trees to gene trees, compare different genes, D statistic determines introgression: ABBA-BABA; if no introgression then D=0 but introgression is there then D is not equal to 0)

     *Limitation:* throws out data, requires many genomes, same value=multiple expectation.

2.   **Liklihood/model based methods**:

                                **Allele frequency spectrum:** <u>is the [distribution](https://en.wikipedia.org/wiki/Frequency_distribution) of the [allele frequencies](https://en.wikipedia.org/wiki/Allele_frequency) of a given set of [loci](https://en.wikipedia.org/wiki/Locus_(genetics)) (often [SNPs](https://en.wikipedia.org/wiki/SNPs)) in a population or sample.

                               Uses count data: distribution with characteristics shape

                               Neutral, bottleneck and selective sweeps.

                                 *Assumptions:*

                                 a. Allele SNPs: independent 

                                 b. Free recombination among SNPs

                                 c. mutation rates are equal

                                 *Limitation:*

                                 a. Loose a lot of data

                                 B. Expensive

                                 **Genealogy sampling:** Multiple regions

                                 *Assumptions-*


     1. Free rcombination among gene
     2. Complete lineage with loci
     3. mutation rates vary across genome
     4. No recombination share common ancestor
    
     **Liklihood free method**: Approx Bayesian comp. (ABC)
    
     Simulations under model of interest (easy)			

#### 3. Historical gene flow + LD patterns

**Distribution of haplotype lengths**

1. recombination leads to shorter fragment over time
2. diificult to ID migrant haplotypes
3. other chemographic events - loci lengths

**Approx of conditional likelihoods**

ancestral recombination graphs (ARG)

*Limitation:*

very complex

difficult to interpret the correct ARG

#### 4. NGS advantages and limitation

**Advantage:**

1. Good estimation  of recombination genes
2. Large area for genomic scans

Limitation:

### Terminal codes

in background

```
[aadas@pbio381 reads2snps]$ screen
```

```
$ cd /data/project_data/snps/reads2snps
$ /data/popgen/dNdSpiNpiS_1.0 -alignment_file=SSW_by24inds.txt.fas -ingroup=sp -out=~/dNdSpiNpiS_output
```

**detach:** Control + A + D

**Reattach** 

```
screen -r
```

While we wait for that to chug along (it'll calculate confidence intervals from 10,000 bootstraps…which takes ~5 hours on our data), we can look at the summary output from the smaller VCF file run previously on just 1 sample library per individual:

```
$ cat SSW_bamlist.txt.sum
```

Selected ingroup species: sp

Number of analyzed individual: 24 (from 1 population(s))

Total number of contig used for sequence analysis: 1113

Total number of SNPs: 5040

- Biallelic: 4991

- Triallelic: 49

- Quadriallelic: 0

##### higher values more homozygous

Fit:

Average Fit: -0.0507419 [-0.06817; -0.031933]

(Fit calculated in 902 contigs)

Weir & Cockerham Fit (Evolution 1984):

Average Weir & Cockerham Fit: 0.00703754 [-0.017669; 0.032047]

(Fit calculated in 902 contigs)

piN/piS ratio:

**Synonymous-**

Average piS in focal species: 0.00585312 [0.005172; 0.006598]

##### Non-synonymous-

*selection* is eliminating bcoz the value is really low

higher values means slelection is not good eliminating the non-syn deletorious mutation

Average piN in focal species: 0.00154546 [0.00133; 0.001782]

Average piN / average piS: 0.264041 [0.223914; 0.310575]

(piS and piN calculated in 902 contigs of average length 50)

##### PiS: 0.00585 and confidence interval [0.005172; 0.006598]

**piN / piS: 0.264041 confidence interval is [0.223914; 0.310575]**

#### Transfer the file

```
#open a new window in your mac
ip0af5268c:users aayudhdas$ cd ~/Dropbox/Aayudh_UVM/ecological\ genomics/pop_gen
# create new folder
ip0af5268c:ecological genomics aayudhdas$ mkdir pop_gen
#now transfer the file from the location, your directory is=pwd
ip0af5268c:pop_gen aayudhdas$ scp aadas@pbio381.uvm.edu:/data/project_data/snps/reads2snps/Romiguier_nature13685-s3.csv .
```

R script

```
setwd("/Users/aayudhdas/Dropbox/Aayudh_UVM/ecological genomics/pop_gen/")
list.files()
# Read in the Romiguier data:
Rom <- read.csv("Romiguier_nature13685-s3.csv", header=T)
# Import OK?
str(Rom) 
head(Rom)
# Looks good!
# Now let's look at how the strength of purifying selection (piN/piS) compares to the size of Ne (piS). We'll plot these on a log scale to linearize the relationship.
plot(log(Rom$piS), log(Rom$piNpiS), pch=21, bg="blue", xlab="log Synonymous Nucleotide Diversity (piS)", ylab="log Ratio of Nonysn to Syn Diversity (piN/piS)", main="Purifying Selection vs. Effective Population Size")
# Now let's add our SSW points to the existing plot and give them a different symbol
points(log(0.00585312), log(0.264041), pch=24, cex=1.5, bg="red") 
# We can also add a regression line to the plot to see how far off the SSW estimates are from expectation
reg <- lm(log(Rom$piNpiS) ~ log(Rom$piS)) # Fits a linear regression
abline(reg) # adds the regression line to the plot
# It would be useful to highlight the other echinoderms in the dataset...do our seastars behave similarly?
echino <- Rom[which(Rom$Phylum=="Echinodermata"),] # subsets the data
points(log(echino$piS), log(echino$piNpiS), pch=21, bg="red") # adds the points
# Lastly, let's add a legend:
legend("bottomleft", cex=1, legend=c("Metazoans", "Echinoderms", "P. ochraceus"), pch=c(21,21,24), col=c("blue", "red", "red"))
# Pisaster seems to be in a group with other echinoderms that have relaxed purifying selection (high piN/piS), given their Ne...Interesting! Can we hypothesize why this might be?
```

## 

<div id='id-section17'/>

### Page 17: 2017-03-27. Population genetics 5

Load libraries:

```
library(vcfR)
library(adegenet)

```

Read the vcf SNP data (the one I unzipped) into R and call it "vcf1"

```
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")

```

I got an error message saying that this file dosn't seem to exist.

Tried again, this time instead of running the version straight up; I copied and pasted my direct file name into the command line (see below):

```
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8..recode.vcf")

```

It worked!

*I edited the script to work with my file

Since we have SOOOO many SNPS we will use one of the packages we loaded (adegenet) to store the large # of SNPS; this lets us analyze the SNPs easier

```
gl1 <- vcfR2genlight(vcf1)
print(gl1)

```

By printing the resulting file you can see that it has the right number of SNPs (5,317) and individuals (24)

We can check the file for certain information using the following commands:

```
gl1$ind.names
gl1$loc.names[1:10]
gl1$chromosome[1:3]

```

The first command (gl1$ind.names) tells you the names of each individual (they just didn't all fit on the same row so they are broken up *# in [] = # individual it starts with

The 2nd command (gl1$loc.names[1:10]) tells you the names of each loci

- Everything to the LEFT of the "_" = transcript ID
- Everything to the RIGHT of the "_" = position (ie 7456=bp location where loci (SNP) occurs within the transcript

The 3rd command (gl1$chromosome[1:3]) shows you the transcript ID

At this point we are adding pop to the field (actually happens at line 46 but we are building up to it)

```
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file

```

The first command read in the meta data and called it "ssw_meta"

The second command organized the meta data so that it was sorted by individual ID

We want to confirm that the IDs are ordered the same in both the SNP data file (the we read in earlier and called "gl1") and the meta data file (that we just read in and called "ssw_meta"

When we match data from ssw_meta to gl1 they just read it in in the order that it is in. That is why we are making sure that the order is the same for each file so that the data will match up correctly

```
gl1$ind.names
ssw_meta$Individual

```

Here is what it looks like after running those commands:

```
> gl1$ind.names
 [1] "03" "07" "08" "09" "10" "14" "15" "19"
 [9] "20" "22" "23" "24" "26" "27" "28" "29"
[17] "31" "32" "33" "34" "35" "36" "37" "38"
> ssw_meta$Individual
 [1]  3  7  8  9 10 14 15 19 20 22 23 24 26 27
[15] 28 29 31 32 33 34 35 36 37 38

```

They don't have to match exactly (07 and 7 is ok in this case) The # in the brackets is only indicating where the line picks up (for example [15] means it picks up with ind #15 (not ind15); it is NOT the row #

Now we will assign locality (first line of code below) and disease status (second line of code below)

- when we say "assign" locality etc we mean we will use tell the gl1 file what the location and disease status of each individual is
- basically we are merging the files to have one file with all the info we need to make the comparisons and calculations we will do later today.

```
gl1$pop <- ssw_meta$Locality 
gl1$other <- as.list(ssw_meta)

```

We can look at the SNP data by generating the following plot:

```
glPlot(gl1, posi="bottomleft")

```

Notes on the plot:

- white space = missing data
- Number of 2nd allele key explained:
- 0= AA
- 1= AT
- 2= TT
- x-axis: each column is a different SNP (there are 5000+ columns)
- y-axis: each row is a different individual

If you want to know which individuals have missing data (white) use the following command to have it tell you what individual is at which row (right now the y axis list the individuals but dosn't correlate with the ind ID (ie row 2 is actually individual 7)

```
gl1$ind.names[2]

```

too look up more then one at a time:

```
gl1$ind.names[c(2, 10)]

```

Now we will calculate THEN plot the PCA (principle coordinate analysis) *NOTE: "pca1" does NOT equal principle componant 1 (its just what we are calling the first pca we are running)

1. Calculate PCA(1) (line 41)

```
pca1 <- glPca(gl1, nf=4, parallel = F)
pca1 

```

nf = # of principle componants to fit in one model second command prints summary

"$" in summary (line 42) show you the lines you can look at (ie scores, loadings etc) if you want to look at specific categories you can run the following command:

```
pca1$scores

```

1. Plot PCA(1) (line 45)

```
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=gl1$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$pop), 
       pch=20, 
       col=c("black", "red"))

```

The above command:

- plots the principle component calcualtions performed by running line 41
  - line 41 calculated PC for ALL SNPs (not just ones specific to health or location, etc)
  - "col=gl1$pop" is the specific part that tells it to color the dots based on location
  - NOTE: by telling it to color them by location it does NOT change where the points are; it mearly colors them based on extra info. The dots are ploted based on PC which was calculated for SNPs and did NOT take these extra variable into acount

Detailed explination of what means what in the command we ran (line 45-53)

```
plot(pca1$scores[,1], pca1$scores[,2], #Calc 4 but only looking at 1 and 2
     cex=2, #point/dot size; will assume cex=1 if you don't specify
	 pch=20 #type of symbol (20=filled circle; can seach for other codes via google or "?pch")
	 col=gl1$pop, #Tell it what to look at (col=color, we said to look at location (which we assigned to "pop")
     xlab="Principal Component 1", #title for x-axis 
     ylab="Principal Component 2", #title for y-axis 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)") #title
legend("topleft", #where to put legend 
       legend=unique(gl1$pop), #if you don't specify "unique" it will list allllll the int and subs (ie int int int sub sub int..); "legend=unique" tells it list each diff variable once (int, sub)
       pch=20, 
       col=c("black", "red")) #telling it to match the colors used in the plot

```

Extra explination of PCA plots (in general):

- the data points in the (+) associate more with the alternative allele
- data points in the (-) associate more with the reference allele (which was determined during the transcriptome assembly (using trinity))
- Each data point represents the AVERAGE of the SNPs per individual (aka its not showing us which SNP associate with which allele; it shows us OVERALL more SNPs associate with the alt/ref allele)

Now we will change the color of the dots to show health staus

- NOTE: remember this will NOT change the data itself but will mearly show you the health status vs location associate with each individual

```
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=gl1$other$Trajectory, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$other$Trajectory), 
       pch=20, 
       col=unique(gl1$other$Trajectory))

```

Now that we have out PCA data we can ask more specific questions like:

*Which SNPs load most strongly on the 1st PC axis?*

```
loadingplot(abs(pca1$loadings[,1]),
            threshold=quantile(abs(pca1$loadings), 0.999))

```

Some notes about the command you just ran above:

- 0.999 = 99.9 percentile (it set the threshold (the horizontal line in the figure) to "cut off" and identify SNPs that were beyond that percentile (?right?)
  - you can try running it at 0.95 (95th percentile) and it will identify MANY more SNPs; this is because it is a less stringent cut off
  - with SNPs since we have so many you need a super stringent cut off (which is why we use 0.999; acceptable = <1%)

If we want to identify (name those SNPs that were above the threshold we can run the following command:

```
gl1$loc.names[which(abs(pca1$loadings)>quantile(abs(pca1$loadings), 0.999))]

```

to look up the actual loading value associated with a particular SNP (for example SNP#3412) we do the following:

```
pca1$loadings[3412]

```

it spat out:

```
[1] 0.1397889

```

**NOW WE MOVE ONTO THE DAPC METHODS**

DAPC (Discriminant analysis of principle componants)

- Tries to decide which SNPs *attributed* to the health statusnot grouping "these SNPs were found in indiv with HH or HS...instead we are trying to find that SNPs that might have caused/that attribute to the health statusgetting more specific then PCA

*The below command does the dapc and will re-run the PCA calc:

```
disease.dapc <- dapc(gl1, pop=gl1$other$Trajectory, n.pca=8, n.da=3,
                     var.loadings=T, pca.info=T, parallel=F)

```

Detailed expliation of each part of the above command:

```
disease.dapc <- dapc(gl1,
				pop=gl1$other$Trajectory, #pop = what you are grouping it by (here we are grouping it by disease)
				n.pca=8, #how many pcs to use; it will calculate 8 here); rule of thumb: don't want to use n/3 pcs (ie for us don't use mroe then # of ind/3; 24/3=8)
				n.da=3, # number of discriminant functions it is using; we have 4 health statuses and you should be able to describe anything with the #variables-1 (ie for us 4-1=3)
                var.loadings=T, # means we want to look at them so output them
				pca.info=T, #we want to beable to ask info abt pc axis in the output 
				parallel=F) #if you need more cores you tell it to use more by making "parallel=T"; here we tell it to use only one because this command can't seem to use more then one 

```

Look at summary of dapc we just ran by running the following:

```
disease.dapc

```

*loadings = info for each pc we calc (should be 8) *PC.loadings = info for all the SNPs

Now we will take the dapc data and plot it; we will use a scatter plot:

```
scatter.dapc(disease.dapc, grp=gl1$other$Trajectory, legend=T)

```

shows all three DA (Discriminant analyses) in the bottom right bar graph style shaded box but it ploted only 1 and 2 (the horizontal and vertical lines that form a "+" in your figure)

We can also show the data as a stacked bar graph

- this helps us to see the probablity for each health status that dacp calculated for each individaul

```
compoplot(disease.dapc)

```

You can use the above plot to compare to the scatter plot generated just before (line 84) to help identify individuals who are linked as (for example) HS but their dot is very close to HH

Finally (for today) you can generate another loading plot to see which SNPs contribute most to distiguishing H vs S individuals

```
loadingplot(abs(disease.dapc$var.load), 
            lab.jitter=1, 
            threshold=quantile(abs(disease.dapc$var.load), probs=0.999))
```

------

<div id='id-section18'/>

### Page 18: 2017-03-29. Population genetics 6

### INFO UPDATE: Detecting local adaptation from population genomic outlier analysis

#### Local adaptation



#### Different approaches

1. Genetic environment association analysis
2. Different outlier method 
   1. FST

#### Common obstcles 

1. Confounding factors
   1. Degrading history
   2. Neutral population structure 
   3. Background selection 
2. Missing genome
   1. Reduced representation 
   2. Missing structural variants in reference
   3. Loss of repetititve regions /paralogs
3. Missing landscape
   1. Low resolution environment data
   2. Scale of local adaptation
   3. Multi colinearity 

#### Solutions

1. Confounding factor
   1. null geographic models
   2. Relatedness among samples
2. Missing genome: exome and RNA seq. WGS and reference genome, depth coverage.
3. Missing landscape: Know systems

#### Other considerations

1. Sampling strategy
2. Multiple comparison
3. Genomic objective

#### Final notes 



### Terminal code

Tutorial 

```

 
/data/project_data/snps/reads2snps/SSW_tidal.pops
/data/project_data/snps/reads2snps/vcf2admixture_SSW.spid
/data/project_data/snps/reads2snps/vcf2geno.sh
The ready-to-go geno file is located on our server here:

 
/data/project_data/snps/reads2snps/SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno
In the same path, you should also see a bash script:

 
/data/project_data/snps/reads2snps/ADMIX.sh
Use cp to copy the .geno and ADMIX.sh files to your home directory on the server, then cd there and confirm the files are present.

From within your home directory, open the ADMIX.sh script in vim. Let's walk through what each step is doing:

 
#!/bin/bash
​
# Run ADMIXTURE to determine the number of genetic clusters in the SNP data, 
# and the ancestry proportions of each individual
​
# Here's the utility of 'for loops'...
​
for K in {1..10}
​
do
​
admixture -C 0.000001 --cv ./SSW_all_biallelic.MAF0.02.Miss1.0.recode.vcf.geno $K \
| tee log${K}.out
​
done
​
# After the for loop finishes, you can use 'grep' to grab the values of the CV from each separate log file and append them into a new summary text file.
​
grep CV log*.out >chooseK.txt
When you're ready to go, exit vim to return to the command line, and execute the script.

 
$ bash ADMIX.sh
```

code

```
  383  2017-03-29 10:33:16 cp ADMIX.sh .
  384  2017-03-29 10:34:38 cp /data/project_data/snps/reads2snps/ADMIX.sh/* . 
  385  2017-03-29 10:35:33 cp -r /data/project_data/snps/reads2snps/* . &
  386  2017-03-29 10:35:56 cd..
  387  2017-03-29 10:36:02 pd
  388  2017-03-29 10:36:04 pwd
  389  2017-03-29 10:36:14 cd ~/
  390  2017-03-29 10:36:28 cp -r /data/project_data/snps/reads2snps/* . &
  391  2017-03-29 10:36:32 ll
  392  2017-03-29 10:38:28 cd /data/project_data/snps/reads2snps/
  393  2017-03-29 10:38:29 ll
  394  2017-03-29 10:38:52 cp SSW_tidal.pops ~/
  395  2017-03-29 10:39:08 cp vcf2admixture_SSW.spid ~/
  396  2017-03-29 10:39:20 cp vcf2geno.sh ~/
  397  2017-03-29 10:39:24 cd ~/
  398  2017-03-29 10:39:27 l
  399  2017-03-29 10:39:28 ll
  400  2017-03-29 10:39:57 gunzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz 
  401  2017-03-29 10:39:59 ll
  402  2017-03-29 10:41:08 vi vcf2geno.sh 
  403  2017-03-29 10:47:08 vi vcf2geno.sh
  404  2017-03-29 10:50:48 ./vcf2geno.sh 
  405  2017-03-29 10:51:05 ll
  406  2017-03-29 10:51:57 vi SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno
  407  2017-03-29 10:54:40 vi vcf2geno.sh
  408  2017-03-29 11:01:26 vi ADMIX.sh 
  409  2017-03-29 11:07:20 bash ADMIX.sh 
  410  2017-03-29 11:07:31 vi ADMIX.sh 
  411  2017-03-29 11:08:13 bash ADMIX.sh 
  412  2017-03-29 11:11:17 cat chooseK.txt 
  413  2017-03-29 11:16:05 vi ADMIX.sh 
  414  2017-03-29 11:16:25 bash ADMIX.sh 
  415  2017-03-29 11:19:02 cat chooseK.txt 
```

------

<div id='id-section19'/>

### Page 19: 2017-03-31. Homework3_population genetics

#### 1. Filtering strategy 1

```
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz
```

After filtering, kept 24 out of 24 Individuals

After filtering, kept 7486938 out of a possible 7486938 Sites

***a) Biallelic vs. multi-allelic SNPs*:** So, SNPs with >2 alleles probably reflect sequence or mapping errors. We also want to get rid of SNPs showing <2 alleles. **b)** ***Minor allele frequency (MAF):*** Sequencing errors are relatively common, but they tend to happen randomly and affect only 1 read at a time. Thus, if we have a SNP that is only seen very rarely, it may be a sequencing error, and should be discarded. For us, the most liberal MAF filters would be 1 allele copy out of the total 2N copies, or 1/48 = 0.02. *Missing data across individuals:* **c)** ***Missing data across individuals:*** Get rid of sites where  fewer than 80% of our samples have data. Missing data is a problem for any analysis, and population genetic statistics can behave oddly (i.e.. become biased) when a lot of individuals are missing data for a given SNP. 

Now if we do 10% Missing data across individuals

```
vcftools --gzvcf SSW_by24inds.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.9 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.9
```

After filtering, kept **3371** out of a possible 7486938 Sites

#### 2. Filtering Strategy 2 (removing individual)

Check the missing individual

```
[aadas@pbio381 homework3]$ vcftools --vcf SSW_all_biallelic.MAF0.02.Miss0.9.recode.vcf --missing-indv --out ~/SSW_class_filtered_missing-indv
```

see the individuals

```
[aadas@pbio381 pipeline2]$ cat SSW_class_filtered_missing0.9-indv.imiss 
INDV	N_DATA	N_GENOTYPES_FILTERED	N_MISS	F_MISS
03	3371	0	7	0.00207654
07	3371	0	752	0.223079
08	3371	0	10	0.00296648
09	3371	0	34	0.010086
10	3371	0	41	0.0121626
14	3371	0	144	0.0427173
15	3371	0	6	0.00177989
19	3371	0	7	0.00207654
20	3371	0	4	0.00118659
22	3371	0	705	0.209137
23	3371	0	46	0.0136458
24	3371	0	421	0.124889
26	3371	0	421	0.124889
27	3371	0	2	0.000593296
28	3371	0	11	0.00326313
29	3371	0	66	0.0195788
31	3371	0	11	0.00326313
32	3371	0	11	0.00326313
33	3371	0	18	0.00533966
34	3371	0	14	0.00415307
35	3371	0	12	0.00355977
36	3371	0	15	0.00444972
37	3371	0	42	0.0124592
38	3371	0	6	0.00177989
```

Individual 07, 14, 22, 24, 26 can be removed.

```
vcftools --gzvcf SSW_by24inds.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.9 --remove-indv 07 --remove-indv 14 --remove-indv 22 --remove-indv 24 --remove-indv 26 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.9.remove07.14.22.24.26
```

After filtering, kept **4762** out of a possible 7486938 Sites

## PCA analysis (R script)

```
library(vcfR)
library(adegenet)
file.choose()
setwd("/Users/aayudhdas/Dropbox/Aayudh_UVM/ecological genomics/homework")
list.files()
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.9.recode.vcf")
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!
# For info, try:
gl1$ind.names
gl1$loc.names[1:10]
gl1$chromosome[1:3]

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file

gl1$ind.names
ssw_meta$Individual

gl1$pop <- ssw_meta$Location # assign locality info

gl1$other <- as.list(ssw_meta) # assign disease status

glPlot(gl1, posi="bottomleft")

pca1 <- glPca(gl1, nf=4, parallel=F) # nf = number of PC axes to retain (here, 4)
pca1 # prints summary
# Plot the individuals in SNP-PCA space, with locality labels:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=gl1$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$pop), 
       pch=20, 
       col=c("black", "red"))
```

------

<div id='id-section20'/>

### Page 20: 2017-04-03. Detecting local adaptation from population genomic outlier analyses 

Concepts

### Inbreeding produces structural popultion 

Fst = (Ht -Ht)/Ht

Fis = [Exp(Hs)-Obs(Hs)]/Exp(Hs)

Fit = (Ht-Hi)/Ht

#### Selective sweeps change allele freq. in pops

#### Empirical p-values created from distilled of previously neutral loci are useful for find Ns

#### Method Out Flank (2015)

Questions:

1. What challenges do outlier detection methods face?
2. How is LD our friend or foe ?

### R script

```
setwd("/Users/aayudhdas/Dropbox/Aayudh_UVM/ecological genomics/karl")
install.packages("devtools")
library("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library("qvalue")
install_github("whitlock/OutFLANK")

library(OutFLANK)
install.packages("vcfR")
library(vcfR)
install.packages("adegenet")
library("adegenet")

ssw.geno_in =read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno", width=rep(1,24))
ssw.geno = t(ssw.geno_in)

ssw_meta = read.table("ssw_healthloc.txt",T)
ssw_meta = ssw_meta[order(ssw_meta$Individual),]
ssw_meta$Trajectory[which(ssw_meta$Trajectory == "MM")] = NA

OF_SNPs = MakeDiploidFSTMat(ssw.geno, locusNames = seq(1, 5317,1), popNames = ssw_meta$Trajectory)
head(OF_SNPs)
OF_out = OutFLANK(FstDataFrame = OF_SNPs, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin = 0.1, NumberOfSamples = 3, qthreshold = 0.1)
OutFLANKResultsPlotter(OF_out, withOutliers = T, NoCorr = T, Hmin = 0.1, binwidth = 0.005, titletext = "Scan for local selection")
outliers = which(OF_out$results$OutlierFlag=="TRUE")
outliers

vcf1 = read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
dim(vcf1)
vcfann = as.data.frame(getFIX(vcf1))
vcfann[outliers,]

```



------

<div id='id-section21'/>

### Page 21: 2017-04-05. Enrichment and annotation

### INFO UPDATE: 

1. After assembly, annotation map reads .sam files using BWA program.
2. Annotation using .fasta file with BLAST program: This will tell us Biological function and match prediction.
3. Annotation: Blast2GO-paid, Brute force, Pipelines-Trinotate

### Terminal code

```
R code in folder
```

how to make table?

| a    | b    | c    |
| ---- | ---- | ---- |
|      |      |      |

Kk

------

<div id='id-section22'/>

### Page 22: 2017-04-10. Metagenomics 1

### Info update

16S ribosomal RNA genes we will be analyzing

Ans the questions like

a. what micobes are present?

b. what are they doing ?

This is called metagenomics

#### Microbiome: 

assemblage of microbial taxa associated with a host or environment.

* pathogens
* Symbiont 
  * mutualist
  * commensals 
* Holobiont 

How to measure micro biome?

* 16s rRNA
* ITS (Internal Transcribed Spacers)

#### SSW projects: Metagenomics 

* genes being expressed by the microbiome

* function

* All RNA or DNA coming from RNA seq or shot-gun DNA seq

  ### Processing: 

  * Sample tissue extract DNA/RNA
  * PCR target gene (16s)
  * Barcode amplified fragment 
  * cluster based on similarity (97%)
  * BLAST OTUs to determine taxonomy.  

#### The big questions

* Is there a "core" micro biome shared in common among many organisms?
* How much diversity is unique to individuals, populations and communities? Heritable?
* Evolution of micro biome—> adaptation of host to environment ?
* Can knowledge of microbiome predict a response to env selection? 


### Terminal code

```
[aadas@pbio381 ~]$ mkdir 16s_analysis
[aadas@pbio381 ~]$ ll
[aadas@pbio381 16s_analysis]$ cd /data/project_data/16s/map.txt
[aadas@pbio381 16s]$ vi map.txt 
[aadas@pbio381 16s_analysis]$ validate_mapping_file.py -m /data/project_data/16s/^Cp.txt -o validate_map -p -b
[aadas@pbio381 16s_analysis]$ validate_mapping_file.py -m /users/a/a/aadas/16s_analysis/map.txt -o validate_map -p -b
[aadas@pbio381 16s_analysis]$ cd validate_map/
[aadas@pbio381 validate_map]$ ll
[aadas@pbio381 16s_analysis]$ multiple_join_paired_ends.py -i /data/project_data/16s/data_files -o ~/16s_analysis/joined --read1_indicator _R1 --read2_indicator _R2
[aadas@pbio381 16s_analysis]$ ll
[aadas@pbio381 16s_analysis]$ cd joined/
[aadas@pbio381 joined]$ ll
[aadas@pbio381 joined]$ bash /data/project_data/16s/remove-underscore.sh
[aadas@pbio381 joined]$ bash /data/project_data/16s/remove-R1.sh
[aadas@pbio381 joined]$ multiple_split_libraries_fastq.py -i ~/16s_analysis/joined -o ~/16s_analysis/filtered -m sampleid_by_file --include_input_dir_path --remove_filepath_in_name  --mapping_indicator /data/project_data/16s/map.txt
```

after 10 min

```
[aadas@pbio381 16s_analysis]$ cd filtered/
[aadas@pbio381 filtered]$ head seqs.fna 
```

Also, run this command on a few samples and make sure the output ‘test’ file contains some sequences

```
extract_seqs_by_sample_id.py -i seqs.fna -o test -s 04-5-05
```

Have a look at the files produces by the pick_open_reference_otus.py. In the interest of time, we won’t go over what all of these files are but the above links gives a great description of what they are and how they were made. We just want to make sure there are the correct number of samples in the biom file and get a sense of the number of reads per sample.

```	
[aadas@pbio381 otu_table]$ biom summarize-table -i /data/project_data/16s/otu_table/otu_table_mc2_w_tax_no_pynast_failures.biom
```

Num samples: 176

Num observations: 93033

Total count: 8362869

Table density (fraction of non-zero values): 0.028

Counts/sample summary:

 Min: 28412.0

 Max: 77866.0

 Median: 47051.500

 Mean: 47516.301

 Std. dev.: 7637.541

 Sample Metadata Categories: None provided

 Observation Metadata Categories: taxonomy

##### phyloseq in R

```
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
```



------

<div id='id-section23'/>

### Page 23: 2017-04-12. Metagenomics 2

### Info update: Host genetic associations with the micro biome

#### Microbiome- 

* digestive
* immune
* vitamins
* xenobiotics
* Pathogens 

#### Heritabilty: 

Proportion of variance in a host trait; Measured across population; 

Explained by genetic (environmental) 

* Directy: Twin pairs
* GAS

#### GWAS

* MB-traits
* SNP data (host)
* 16s rRNA (MB)

#### Cross study comparisons

* Mice QTL
  * immune flexibility
* Plants QTL
  * species niches and abundances of individual taxa
* Flies 
  * GWAS
  * Barrier defense (digestive)
  * Immunity

#### Limitations 

1. Small sample size so not significant result 
2. Lacking taxa
3. Meta analysis results are inconclusive

### Terminal code

##### 1. Identify chimeric OTUs

```
[aadas@pbio381 ~]$ cd 16s_analysis/
[aadas@pbio381 16s_analysis]$ vsearch --uchime_ref /data/project_data/16s/otu_table/rep_set.fna --chimeras ~/16s_analysis/mc2_w_tax_no_pynast_failures_chimeras.fasta --db /usr/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta 
```

##### 2. Remove these chimeric OTUs from OTU table

```
[aadas@pbio381 16s_analysis]$ filter_otus_from_otu_table.py -i /data/project_data/16s/otu_table/otu_table_mc2_w_tax_no_pynast_failures.biom -o otu_table_mc2_w_tax_no_pynast_failures_no_chimeras.biom -e ~/16s_analysis/mc2_w_tax_no_pynast_failures_chimeras.fasta
```

##### 3. Remove these chimeric OTUs from rep_set_aligned and re-make the phylogenetic tree

```
[aadas@pbio381 16s_analysis]$ filter_fasta.py -f /data/project_data/16s/otu_table/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta -o ~/16s_analysis/rep_set_aligned_pfiltered_no_chimeras.fasta -a ~/16s_analysis/mc2_w_tax_no_pynast_failures_chimeras.fasta -n
```

```
[aadas@pbio381 16s_analysis]$ make_phylogeny.py -i ~/16s_analysis/rep_set_aligned_pfiltered_no_chimeras.fasta -o ~/16s_analysis/rep_set_no_chimeras.tre
^C[aadas@pbio381 16s_analysis]$ screen
[detached from 45250.pts-6.pbio381]
```

#### Frequency filtering (##How many OTUs are left?)

Now we will get rid of low frequency OTUs. We will get rid of any OTU with fewer than 50 total counts across all samples (–min_count 50) as well as any OTU that is in fewer than 25% of samples (–min_samples 44)

```
[aadas@pbio381 16s_analysis]$ filter_otus_from_otu_table.py -i otu_table_mc2_w_tax_no_pynast_failures_no_chimeras.biom -o otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom --min_count 50 --min_samples 44

##How many OTUs are left?
[aadas@pbio381 16s_analysis]$ biom summarize-table -i otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom
```

Num samples: 176

Num observations: 1064

Total count: 7221141

Table density (fraction of non-zero values): 0.538

Counts/sample summary:

 Min: 18858.0

 Max: 66754.0

 Median: 40720.500

 Mean: 41029.210

 Std. dev.: 8381.006

 Sample Metadata Categories: None provided

 Observation Metadata Categories: taxonomy

#### Core diversity analyses in QIIME

```
core_diversity_analyses.py -o core_diversity_filtered -i otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom -m ~/Po/MiSeq/joined/map.txt -t rep_set_no_chimeras.tre -e 20000 -a 8
```

This command also takes a long time so we’ll just look at the output online: [http://www.uvm.edu/~mlloyd/mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered_core_diversity/](http://www.uvm.edu/~mlloyd/mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered_core_diversity/)



------

<div id='id-section24'/>

### Page 24: 2017-04-17. Metagenomics 3

Fred's talk

------

<div id='id-section25'/>

### Page 25: 2017-04-19. Metagenomics 4

### Info update: Raerefying Microbiome Data in NOT acceptable

#### Rarefying and Rarefaction

Past method

1. Proportion
2. High rate of false positives
3. Heteroscedasticity
4. Reduced power

##### Rarefying 

-steps

1. Small N and Nmin
2. Discard Libs with fewer than Lmin
3. Sub-sample longer Ns without replacement

What this does?

* Normaility data
* Removing individual no. of reads

Problems with Raerifying 

* High rate of false positives
* Requires emission actual data
* Reduces statistical power

#### Mixture model

* Not decreasing power 
* Increasing accuracy over other method 

#### Distribution

Binomial with o-inflated gaussian 

accounts for biological variability

* edgeR, DESeq2, phyloseq

### Terminal code

PICRUST can only work with OTUS that were picked by closed reference OTU picking.

How many closed ref OTUs?

```
[aadas@pbio381 ~]$ cd 16s_analysis/
[aadas@pbio381 16s_analysis]$ filter_otus_from_otu_table.py -i ~/16s_analysis/otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom -o ~/16s_analysis/closed_otu_table.biom --negate_ids_to_exclude -e /usr/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta 
[aadas@pbio381 16s_analysis]$ biom summarize-table -i closed_otu_table.biom
```

Num samples: 176
Num observations: **259**
Total count: 1493357
Table density (fraction of non-zero values): 0.546

Counts/sample summary:
 Min: 152.0
 Max: 32426.0
 Median: 6333.500
 Mean: 8484.983
 Std. dev.: 6937.648
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

normalize by picrust’s method

```
normalize_by_copy_number.py -i ~/16s_analysis/closed_otu_table.biom -o ~/16s_analysis/closed_otu_table_norm.biom
```

The output of PICRUST is the same format as an OTU table but instead of OTUS, KEGG Orthology terms are used. For a tab delimited output:

```
predict_metagenomes.py -f -i ~/16s_analysis/closed_otu_table_norm.biom -o ~/16s_analysis/metagenome_predictions.txt -a nsti_per_sample.txt 
```

For a .biom output:

```
predict_metagenomes.py -i ~/16s_analysis/closed_otu_table_norm.biom -o ~/16s_analysis/metagenome_predictions.biom -a nsti_per_sample.txt
```

We can take a look at the tab delimited otu table

```
head metagenome_predictions.txt
```

And count the rows to see how many KO were predicted

```
wc -l metagenome_predictions.txt
```

**6910 metagenome_predictions.txt**

What do you think about the number of KOs predicted considering the number of OTUs we’re using in this analysis? 

Collapse to higher KO hierarchy term

```
categorize_by_function.py -f -i metagenome_predictions.biom -c KEGG_Pathways -l 3 -o metagenome_predictions.L3.txt
```

We also want this output in .biom format to use in R

```
categorize_by_function.py -i metagenome_predictions.biom -c KEGG_Pathways -l 3 -o metagenome_predictions.L3.biom
```

Files I have

-rw-r--r--.  1 aadas users 230K Apr 19 10:42 closed_otu_table.biom

-rw-r--r--.  1 aadas users 254K Apr 19 10:50 closed_otu_table_norm.biom

-rw-r--r--.  1 aadas users 6.2M Apr 19 10:52 metagenome_predictions.txt

-rw-r--r--.  1 aadas users 5.5M Apr 19 10:56 metagenome_predictions.biom

-rw-r--r--.  1 aadas users 6.6K Apr 19 10:56 nsti_per_sample.txt

-rw-r--r--.  1 aadas users 397K Apr 19 11:01 metagenome_predictions.L3.txt

-rw-r--r--.  1 aadas users 417K Apr 19 11:01 metagenome_predictions.L3.biom

```
setwd("~/Dropbox/Aayudh_UVM/ecological genomics/Metagenomics")
library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("ggplot2")
theme_set(theme_bw())
install.packages("RJSONIO")
file.choose()
install.packages("/Users/aayudhdas/Downloads/biom_0.3.12.tar.gz", repos = NULL, type = "source")
library("biom")
setwd("~/Dropbox/Aayudh_UVM/ecological genomics/Metagenomics")
x = read_biom("metagenome_predictions.L3.biom")
otumat = as(biom_data(x), "matrix")
OTU = otu_table(otumat, taxa_are_rows=TRUE)


mapping <- import_qiime_sample_data(mapfilename = 'R_map.txt')

phylo <- merge_phyloseq(OTU, mapping)
phylo

###############################################################################
###Test for DE KO terms between individuals that got sick and those that didn't
###############################################################################

final_pheno = phyloseq_to_deseq2(phylo, ~ Final_phenotype)
final_pheno_results = DESeq(final_pheno, test="Wald")
final_pheno_res = results(final_pheno_results)
summary(final_pheno_res)
head(final_pheno_res)

alpha = 0.05
final_pheno_sigtab = final_pheno_res[which(final_pheno_res$padj < alpha), ]
final_pheno_sigtab= cbind(as(final_pheno_sigtab, "data.frame"), as(tax_table(phylo)[rownames(final_pheno_sigtab), ], "matrix"))
head(final_pheno_sigtab)
final_pheno_sigtab
write.table(final_pheno_sigtab, "Final_pheno_L3.txt", sep="\t")
```



------

<div id='id-section26'/>

### Page 26: 2017-04-24. Project code

Transcoder is-

/data/popgen/TransDecoder-3.0.1

**So, how many genes related to immune response?**

[aadas@pbio381 ~]$ grep "immune" poch_uniprot_GO_nr.txt | wc -l

279

[aadas@pbio381 ~]$ grep "immun*" poch_uniprot_GO_nr.txt | wc -l

321

**put it in a text file**

grep "immun*" poch_uniprot_GO_nr.txt > poch_uniprot_GO_nr_immune.txt

**Now you need the header**

head poch_uniprot_GO_nr.txt -n 1 > head_poch_uniprot_GO_nr_immune.txt

Make a file including the header

cat head_poch_uniprot_GO_nr_immune.txt poch_uniprot_GO_nr_immune.txt > immune.txt

**Now from UPS and ARMS paper FASTA file**

[aadas@pbio381 ~]$ /data/popgen/TransDecoder-3.0.1/TransDecoder.LongOrfs -t Phel_transcriptome_clc_v3.fasta

**make the database** 

```
makeblastdb -in Pycnopodia_longest_orfs.pep -dbtype prot
```

**Script to work on**

```
#!/bin/bash

blastp -query ~/Pisaster_longest_orfs.pep \
       -db ~/Pycnopodia_longest_orfs.pep \
       -out ~/Pisaster_vs_Pycnopodia.outfmt6 \
       -outfmt 6 \
       -evalue 1e-3 \
       -max_target_seqs 1
```

[aadas@pbio381 ~]$ wc -l Pisaster_vs_Pycnopodia.outfmt6

12916 Pisaster_vs_Pycnopodia.outfmt6

R code

```
setwd("~/Dropbox/Aayudh_UVM/ecological genomics/project")
read.table("poch_uniprot_GO_nr.txt")
read.table("Pisaster_vs_Pycnopodia.txt")
PivsPyc <- read.delim('Pisaster_vs_Pycnopodia.txt', header=FALSE)
head(PivsPyc)
Pisaster <- read.delim('poch_uniprot_GO_nr.txt', header=TRUE)
head(Pisaster)
```

Box plot code:

```
setwd("~/Dropbox/Aayudh_UVM/ecological genomics/project")
library("ggplot2")
stressgenes <- read.delim("stress.final.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)
countData <- data.frame(gene=row.names(countData),countData)
head(countData)
str(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
conds$individual<-as.factor(rownames(conds))
head(conds)
conds$day<-as.numeric(substr(conds$day,4,5))

colData <- as.data.frame(conds)
head(colData)

dim(colData)
dim(conds)

library(tidyr)
data_long <- gather(countData, individual, counts, I03_5.08_S_2:I38_6.12_H_0, factor_key=TRUE)
head(data_long)
dim(data_long)
#install.packages("dplyr")
library(dplyr)
names(conds)
alldata<-inner_join(data_long,conds,by="individual")
names(alldata)[1]<- "trinID"
dim(alldata)
head(alldata)
head(stressgenes)
stress_final<-inner_join(alldata,stressgenes,by="trinID")
dim(stress_final)


ggplot(stress_final,aes(x=day,y=log(counts+1),colour=health))+geom_jitter()
ggplot(stress_final,aes(x=health,y=log(counts+1),colour=individual))+geom_jitter()

ggplot(stress_final,aes(x=factor(day),y=log(counts+1),colour=health))+geom_boxplot()
#+geom_line()
ggplot(stress_final,aes(x=factor(health),y=log(counts+1),colour=score))+geom_boxplot()
plotPCA(stress_final)
```

Differential gene expression analysis (Team Sherlock)

```

library(DESeq2)

#Subsetting data
conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1, sep="")
colData1 <- as.data.frame(conds)
colDataDay6<-subset(colData1, colData1$day=="day06")
colDataDay6
colData<-subset(colDataDay6, colDataDay6$indiv != "37")
colData # 7 sick # 11 healthy
# Day 15 H: 10, 24, 27, 31, 33, 34
# Day 15 S: 8, 9, 15, 19, 20, 23
# no 7, 22 (sick) in day 6
#sample log2foldchange of 279 genes with replacement - calculate average 10000 times to generate distribution; which proportion are sig on their own

colDataDay15<-subset(colData1, colData1$day=="day15")
colDataDay15
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
countData<-countData[, which(colnames(countData) %in% row.names(colData))]

#DESeq Model: Day 6, no MMs, 18 individuals (7 sick, 11 healthy); 
library(DESeq2)
ddsFULL <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)

ddsFULL <- ddsFULL[ rowSums(counts(ddsFULL)) > 100, ]

colData(ddsFULL)$health <- factor(colData(ddsFULL)$health, levels=c("H","S"))

ddsFULL <- DESeq(ddsFULL)
resFULL <- results(ddsFULL)
resFULL <- resFULL[order(resFULL$padj),] #sorts according to pvalue
head(resFULL)
summary(resFULL)
names(resFULL)

#immune related genes (269/final 254)
immuneNames <- read.delim("immune.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(immuneNames) <- immuneNames$trinID
immuneRes<-resFULL[which(row.names(resFULL) %in% row.names(immuneNames)),]
dim(immuneRes) # 254 matched -  less than 321 because of low read counts

#Reading in stress related genes (276/ final 271)
stressNames <- read.delim("Pisaster_stress.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(stressNames) <- stressNames$trinID
stressRes<-resFULL[which(row.names(resFULL) %in% row.names(stressNames)),]
dim(stressRes) 

#Reading in virus related genes (325/final 323)
virusNames <- read.delim("virus.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(virusNames) <- virusNames$trinID
virusRes<-resFULL[which(row.names(resFULL) %in% row.names(virusNames)),]
dim(virusRes) 

# Reading in bacteria related genes (1320/final 1299)
bacteriaNames <- read.delim("bacteria.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(bacteriaNames) <- bacteriaNames$trinID
bacteriaRes<-resFULL[which(row.names(resFULL) %in% row.names(bacteriaNames)),]
dim(bacteriaRes)

# Random distribution of full data set 
oPar <- par(no.readonly=TRUE)

par(mfrow=c(1,4))
par(mar=c(1,1,1,1))

# Immune Random null distribution
logfoldrandoMean0 <- vector(mode="numeric")
for(i in 1:10000){
    logfoldrando0<-vector(mode='numeric')
    logfoldrando0 <- sample(x=resFULL$log2FoldChange, size=dim(immuneRes)[1], replace=T)
    logfoldrandoMean0[i]<-mean(logfoldrando0)
}

#Histogram immune
hist(logfoldrandoMean0, 
     breaks=50, 
     main="IMMUNE GENES",
     xlab=NULL,
     ylim=c(0,800),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
)

Interval0 <- quantile(x=logfoldrandoMean0, probs=c(0.025,0.975))
abline(v = Interval0,col="black",lwd=2,lty="dotted")
abline(v=mean(immuneRes$log2FoldChange), col="red",lwd=2)
t.test(x=logfoldrandoMean0, y=immuneRes$log2FoldChange) #pvalue 0.1403

#Stress Random null distribution
logfoldrandoMean1 <- vector(mode="numeric")
for(i in 1:10000){
    logfoldrando1<-vector(mode='numeric')
    logfoldrando1 <- sample(x=resFULL$log2FoldChange, size=dim(stressRes)[1], replace=T)
    logfoldrandoMean1[i]<-mean(logfoldrando1)
}

#Histogram - stress
hist(logfoldrandoMean1, 
     breaks=70, 
     main="STRESS GENES",
     xlab=NULL,
     ylim=c(0,800),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
)

Interval1 <- quantile(x=logfoldrandoMean1, probs=c(0.025,0.975))
abline(v = Interval1,col="black",lwd=2,lty="dotted")
abline(v=mean(stressRes$log2FoldChange), col="red",lwd=2)
t.test(x=logfoldrandoMean1, y=stressRes$log2FoldChange, alternative="two.sided") #p-value = 0.05001

# Stress Random null distribution
logfoldrandoMean2 <- vector(mode="numeric")
for(i in 1:10000){
    logfoldrando2<-vector(mode='numeric')
    logfoldrando2 <- sample(x=resFULL$log2FoldChange, size=dim(virusRes)[1], replace=T)
    logfoldrandoMean2[i]<-mean(logfoldrando2)
}

#Histogram - virus
hist(logfoldrandoMean2,
     breaks=50,
     main="VIRUS GENES",
     xlab="log(fold change)",
     ylim=c(0,900),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
)

Interval2 <- quantile(x=logfoldrandoMean2, probs=c(0.025,0.975))
abline(v = Interval2,col="black",lwd=2,lty="dotted")
abline(v=mean(virusRes$log2FoldChange), col="red",lwd=2)
t.test(x=logfoldrandoMean2, y=virusRes$log2FoldChange) #p-value = 0.5967

# Stress Random null distribution
logfoldrandoMean3 <- vector(mode="numeric")
for(i in 1:10000){
    logfoldrando3<-vector(mode='numeric')
    logfoldrando3 <- sample(x=resFULL$log2FoldChange, size=dim(bacteriaRes)[1], replace=T)
    logfoldrandoMean3[i]<-mean(logfoldrando3)
}

#Histogram - bacteria
hist(logfoldrandoMean3, 
     breaks=50,
     main="BACTERIA GENES",
     xlab="log(fold change)",
     ylim=c(0,900),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2) 


Interval3 <- quantile(x=logfoldrandoMean3, probs=c(0.025,0.975))
abline(v = Interval3,col="black",lwd=2,lty="dotted")
abline(v=mean(bacteriaRes$log2FoldChange), col="red",lwd=2)
t.test(x=logfoldrandoMean3, y=bacteriaRes$log2FoldChange, alternative="two.sided") #p-value = 0.03746

pvalue <- function(nullDat=logfoldrandoMean1, obsDat=stressRes){
    
    obsDat <- obsDat[,2]
    
    xbar <- mean(obsDat)
    mu <- mean(nullDat)
    n <- length(obsDat)
    SD <- sd(obsDat)
    SE <- SD/sqrt(n)
    
    z <- (xbar - mu)/(SE)
    
    pval <- 2*pnorm(-abs(z), lower.tail = TRUE)
    
    out <- list(xbar=xbar, mu=mu, n=n, SD=SD, SE=SE, z=z, pval=pval)
    
    return(out)
}

running test (pvalue) for each of 6 groups


immune <- pvalue(nullDat=logfoldrandoMean0, obsDat=immuneRes)
stress <- pvalue(nullDat=logfoldrandoMean1, obsDat=stressRes)
virus <- pvalue(nullDat=logfoldrandoMean2, obsDat=virusRes)
bacteria <- pvalue(nullDat=logfoldrandoMean3, obsDat=bacteriaRes)

# creating data frame for figure: (genegroup, mean, SE etc.)
GeneGroup <- c("immune",
               "stress",
               "virus",
               "bacteria")


xBar <- c(immune$xbar,
          stress$xbar,
          virus$xbar,
          bacteria$xbar)


SE <- c(immune$SE,
        stress$SE,
        virus$SE,
        bacteria$SE)

pVal <- c(immune$pval,
          stress$pval,
          virus$pval,
          bacteria$pval)


DF3 <- data.frame(GeneGroup, xBar, SE, pVal)

# calculating confidence interval from null data:
meanCI <- mean(immune$mu,
               stress$mu,
               virus$mu,
               bacteria$mu)

# caclulating quantiles foe each group
IntervalIM <- quantile(x=logfoldrandoMean0, probs=c(0.025,0.975))
IntervalST<- quantile(x=logfoldrandoMean1, probs=c(0.025,0.975))
IntervalVI <- quantile(x=logfoldrandoMean2, probs=c(0.025,0.975))
IntervalBA <- quantile(x=logfoldrandoMean3, probs=c(0.025,0.975))

lower <- mean(IntervalIM[1],IntervalST[1],IntervalVI[1],IntervalBA[1])
upper <- mean(IntervalIM[2],IntervalST[2],IntervalVI[2],IntervalBA[2])

# plotting using ggplot

#choosing color pallet
colors <- c(rep("dodgerblue4",4))

#Create a bar graph for with CI and SE bars
plot1 <- ggplot(DF3, aes(x=GeneGroup,
                         y=xBar,
                         fill=GeneGroup)) + 
    geom_bar(stat="identity",
             color = "black",
             position=position_dodge()) + labs(x="Gene Group", y = "log(Fold Change)") + geom_errorbar(aes(ymin=xBar-SE, ymax=xBar+SE), width=.2, position=position_dodge(.9)) 

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors) + coord_cartesian(ylim = c(-.3, .3)) + annotate("segment", x = .5,xend = 4.5 ,y = meanCI, yend = meanCI, col = "red") + annotate("segment", x = .5, xend = 4.5, y = upper, yend = upper, col = "red", linetype = 2, size = 0.3) + annotate("segment", x = .5, xend = 4.5, y = lower, yend = lower, col = "red", linetype = 2, size = 0.3) + annotate("segment", x = .5, xend = 4.5, y = 0, yend = 0, col = "black", size = 1) + annotate(geom = "text", x = 3, y = .2, label = "*",cex = 8) + theme(legend.position=c(3, 3))


# Normalizing count data
library(DESeq2)
dds2 <- estimateSizeFactors(ddsFULL)
normcounts<-as.data.frame(counts(dds2, normalized=TRUE))

# Immune-related counts
immuneNormCounts<-normcounts[which(row.names(normcounts) %in% row.names(immuneNames)),]
immuneNorm<-as.data.frame(immuneNormCounts)
hstat<-c("S","S","H","S","H","H","H","H","S","H","S","S","H","H","H","H","H","S")
library(data.table)
immuneNormList<-as.list(immuneNorm)
immNorm<-rbindlist(list(immuneNormList))
timmNorm<-as.data.frame(t(immNorm))

## Random Normalized Counts
randomNormCounts<-normcounts[sample(row.names(normcounts), 254),]
randNorm<-as.data.frame(randomNormCounts)
hstat<-c("S","S","H","S","H","H","H","H","S","H","S","S","H","H","H","H","H","S")
library(data.table)
randList<-as.list(randNorm)
randNorm2<-rbindlist(list(randList))
tRandNorm<-as.data.frame(t(randNorm2))

random<-as.data.frame(tRandNorm)
randomMeans<-apply(random, 1, mean)

immune<-as.data.frame(timmNorm)
immuneMeans<-apply(immune,1,mean)

genes<-rep(c("Immune","Random"), each=18)
means<-as.data.frame(c(immuneMeans,randomMeans))
means$genes<-genes
names(means)[1]<-"counts"
summary(aov(counts~genes, data=means))
par(mar=c(5,5,5,5))
boxplot(counts~genes, data=means, col=c("goldenrod","forestgreen"),ylab="Normalized read counts", cex=1, cex.axis=2, cex.lab=2)


#Volcano Plot
par(mfrow=c(1,1))
library(DESeq2)
with(resFULL, plot(log2FoldChange, -log10(padj), pch=20, main="", xlim=c(-6,5),cex.main=1.8, cex.axis=1.5, cex.lab=1.5))
with(subset(resFULL, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(resFULL, padj<.1 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(resFULL, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))
max(resFULL$log2FoldChange)
library(calibrate)
resFULL$Gene<-substr(row.names(resFULL), 9,21) #Made a new column so I can label genes
with(subset(resFULL, padj<.005 & abs(log2FoldChange)>2), textxy(log2FoldChange, -log10(padj), labs=Gene, cex=0.8))
```

