### Page 19: 2017-03-31. Homework3_population genetics

#### 1. Filtering

```
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz
```

After filtering, kept 24 out of 24 Individuals

After filtering, kept 7486938 out of a possible 7486938 Sites

***a) Biallelic vs. multi-allelic SNPs*:** So, SNPs with >2 alleles probably reflect sequence or mapping errors. We also want to get rid of SNPs showing <2 alleles. **b)** ***Minor allele frequency (MAF):*** Sequencing errors are relatively common, but they tend to happen randomly and affect only 1 read at a time. Thus, if we have a SNP that is only seen very rarely, it may be a sequencing error, and should be discarded. For us, the most liberal MAF filters would be 1 allele copy out of the total 2N copies, or 1/48 = 0.02. *Missing data across individuals:* **c)** ***Missing data across individuals:*** Get rid of sites where  fewer than 80% of our samples have data. Missing data is a problem for any analysis, and population genetic statistics can behave oddly (i.e.. become biased) when a lot of individuals are missing data for a given SNP. 

```
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/homework3 SSW_all_biallelic.MAF0.02.Miss0.8
```

After filtering, kept 5317 out of a possible 7486938 Sites

**d) Hardy-Weinberg equilibrium expectations**: **(1=p^2 + 2pq + q^2)**. Use the quality-filtered file we generated above as input.

```
[aadas@pbio381 homework3]$ gzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf 
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --hardy
```

After filtering, kept 24 out of 24 Individuals

Outputting HWE statistics (but only for biallelic loci)

​	HWE: Only using fully diploid SNPs.

After filtering, kept 5317 out of a possible 5317 Sites

Now after 4 filters lets check in R

```
> hwe<-read.table("out.hwe", header=T)
> str(hwe)
'data.frame':	1494 obs. of  8 variables:
 $ CHR               : Factor w/ 360 levels "TRINITY_DN27892_c0_g1_TRINITY_DN27892_c0_g1_i1_g.3123_m.3123",..: 355 355 355 355 355 355 355 355 355 355 ...
 $ POS               : int  4733 5850 5865 5869 5874 6096 6146 6201 6289 6325 ...
 $ OBS.HOM1.HET.HOM2.: Factor w/ 53 levels "10/10/4","10/11/3",..: 43 52 43 43 43 43 43 36 43 53 ...
 $ E.HOM1.HET.HOM2.  : Factor w/ 22 levels "10.01/10.98/3.01",..: 17 2 17 17 17 17 17 14 17 3 ...
 $ ChiSq_HWE         : num  0.0109 2.3438 0.0109 0.0109 0.0109 ...
 $ P_HWE             : num  1 0.2 1 1 1 ...
 $ P_HET_DEFICIT     : num  1 0.979 1 1 1 ...
 $ P_HET_EXCESS      : num  1 0.164 1 1 1 ...
```

```
>summary(hwe)
TRINITY_DN45147_c0_g1_TRINITY_DN45147_c0_g1_i3_g.18680_m.18680:  32  
 TRINITY_DN46382_c0_g1_TRINITY_DN46382_c0_g1_i1_g.22149_m.22149:  29  
 TRINITY_DN45750_c0_g1_TRINITY_DN45750_c0_g1_i2_g.20209_m.20209:  26  
 TRINITY_DN47302_c3_g1_TRINITY_DN47302_c3_g1_i2_g.25471_m.25471:  23  
 TRINITY_DN45186_c3_g1_TRINITY_DN45186_c3_g1_i3_g.18787_m.18787:  18  
 TRINITY_DN46789_c1_g3_TRINITY_DN46789_c1_g3_i1_g.23393_m.23393:  18  
 (Other)                                                       :1348  
      POS         OBS.HOM1.HET.HOM2.        E.HOM1.HET.HOM2.
 Min.   :   1.0   23/1/0 :849        23.01/0.98/0.01:849    
 1st Qu.: 177.0   22/2/0 :192        22.04/1.92/0.04:231    
 Median : 318.0   21/3/0 :106        21.09/2.81/0.09:114    
 Mean   : 614.7   20/4/0 : 81        20.17/3.67/0.17: 89    
 3rd Qu.: 715.5   19/5/0 : 49        19.26/4.48/0.26: 60    
 Max.   :6511.0   23/0/1 : 39        18.38/5.25/0.38: 40    
                  (Other):178        (Other)        :111    
   ChiSq_HWE             P_HWE          P_HET_DEFICIT        P_HET_EXCESS      
 Min.   : 0.000086   Min.   :1.00e-07   Min.   :0.0000001   Min.   :0.0006658  
 1st Qu.: 0.010865   1st Qu.:1.00e+00   1st Qu.:1.0000000   1st Qu.:0.9787234  
 Median : 0.010865   Median :1.00e+00   Median :1.0000000   Median :1.0000000  
 Mean   : 0.995220   Mean   :9.26e-01   Mean   :0.9288662   Mean   :0.9469668  
 3rd Qu.: 0.106667   3rd Qu.:1.00e+00   3rd Qu.:1.0000000   3rd Qu.:1.0000000  
 Max.   :24.000000   Max.   :1.00e+00   Max.   :1.0000000   Max.   :1.0000000  
```

```
> which(hwe$P_HET_DEFICIT<0.01)
[1] 1014 1034 1036 1321 1323 1344 1438 1440
> hwe[which(hwe$P_HET_DEFICIT<0.01),]
                                                                  CHR POS
1014 TRINITY_DN45155_c27_g2_TRINITY_DN45155_c27_g2_i2_g.18743_m.18743 216
1034 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742  99
1036 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742 138
1321     TRINITY_DN39079_c3_g1_TRINITY_DN39079_c3_g1_i1_g.8354_m.8354 244
1323     TRINITY_DN39079_c3_g1_TRINITY_DN39079_c3_g1_i1_g.8354_m.8354 279
1344     TRINITY_DN39696_c4_g1_TRINITY_DN39696_c4_g1_i1_g.8926_m.8926 283
1438   TRINITY_DN42225_c1_g1_TRINITY_DN42225_c1_g1_i1_g.12458_m.12458 220
1440   TRINITY_DN42225_c1_g1_TRINITY_DN42225_c1_g1_i1_g.12458_m.12458 255
     OBS.HOM1.HET.HOM2. E.HOM1.HET.HOM2. ChiSq_HWE        P_HWE P_HET_DEFICIT
1014             22/0/2  20.17/3.67/0.17        24 1.418440e-03  1.418440e-03
1034            11/0/13  5.04/11.92/7.04        24 9.114786e-08  9.114786e-08
1036             19/0/5  15.04/7.92/1.04        24 6.498371e-06  6.498371e-06
1321            13/0/11  7.04/11.92/5.04        24 9.114786e-08  9.114786e-08
1323             22/0/2  20.17/3.67/0.17        24 1.418440e-03  1.418440e-03
1344            13/0/11  7.04/11.92/5.04        24 9.114786e-08  9.114786e-08
1438            11/0/13  5.04/11.92/7.04        24 9.114786e-08  9.114786e-08
1440             22/0/2  20.17/3.67/0.17        24 1.418440e-03  1.418440e-03
     P_HET_EXCESS
1014            1
1034            1
1036            1
1321            1
1323            1
1344            1
1438            1
1440            1
```

#### p value

```
> hwe[which(hwe$P_HET_EXCESS<0.001),]
                                                               CHR POS
588 TRINITY_DN44774_c5_g1_TRINITY_DN44774_c5_g1_i1_g.17737_m.17737 138
    OBS.HOM1.HET.HOM2. E.HOM1.HET.HOM2. ChiSq_HWE        P_HWE P_HET_DEFICIT
588             4/20/0  8.17/11.67/4.17   12.2449 0.0006987037             1
    P_HET_EXCESS
588 0.0006657733
```

```
> hwe[which(hwe$P_HET_DEFICIT<0.001),]
                                                                  CHR POS
1034 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742  99
1036 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742 138
1321     TRINITY_DN39079_c3_g1_TRINITY_DN39079_c3_g1_i1_g.8354_m.8354 244
1344     TRINITY_DN39696_c4_g1_TRINITY_DN39696_c4_g1_i1_g.8926_m.8926 283
1438   TRINITY_DN42225_c1_g1_TRINITY_DN42225_c1_g1_i1_g.12458_m.12458 220
     OBS.HOM1.HET.HOM2. E.HOM1.HET.HOM2. ChiSq_HWE        P_HWE P_HET_DEFICIT
1034            11/0/13  5.04/11.92/7.04        24 9.114786e-08  9.114786e-08
1036             19/0/5  15.04/7.92/1.04        24 6.498371e-06  6.498371e-06
1321            13/0/11  7.04/11.92/5.04        24 9.114786e-08  9.114786e-08
1344            13/0/11  7.04/11.92/5.04        24 9.114786e-08  9.114786e-08
1438            11/0/13  5.04/11.92/7.04        24 9.114786e-08  9.114786e-08
     P_HET_EXCESS
1034            1
1036            1
1321            1
1344            1
1438            1
```



```
[aadas@pbio381 homework3]$ vcftools --vcf biallelic.MAF0.2.Miss0.8.recode.vcf --geno-r2
[aadas@pbio381 homework3]$ vim out.geno.ld
> LD <- read.table("out.geno.ld", header=T)
> str(LD)
'data.frame':	29611 obs. of  5 variables:
 $ CHR   : Factor w/ 737 levels "TRINITY_DN29134_c0_g1_TRINITY_DN29134_c0_g1_i1_g.3467_m.3467",..: 724 724 724 724 724 724 724 724 724 724 ...
 $ POS1  : int  5715 5715 5715 5715 5715 5715 5715 5715 5715 5715 ...
 $ POS2  : int  5723 5726 5728 6313 6354 6396 6481 6647 6731 6780 ...
 $ N_INDV: int  21 21 21 19 18 19 21 19 20 20 ...
 $ R.2   : num  0.00526 0.00526 0.00526 0.02206 0.00735 ...
  $ R.2   : num  0.00526 0.00526 0.00526 0.02206 0.00735 ...
> LD$dist<- abs(LD$POS1-LD$POS2)
> str(LD)
'data.frame':	29611 obs. of  6 variables:
 $ CHR   : Factor w/ 737 levels "TRINITY_DN29134_c0_g1_TRINITY_DN29134_c0_g1_i1_g.3467_m.3467",..: 724 724 724 724 724 724 724 724 724 724 ...
 $ POS1  : int  5715 5715 5715 5715 5715 5715 5715 5715 5715 5715 ...
 $ POS2  : int  5723 5726 5728 6313 6354 6396 6481 6647 6731 6780 ...
 $ N_INDV: int  21 21 21 19 18 19 21 19 20 20 ...
 $ R.2   : num  0.00526 0.00526 0.00526 0.02206 0.00735 ...
 $ dist  : int  8 11 13 598 639 681 766 932 1016 1065 ...
> pdf("LD.plot.pdf")
> plot(LD$dist,LD$R.2)
> dev.off()
```

