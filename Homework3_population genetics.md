### Page 19: 2017-03-31. Homework3_population genetics

#### 1. Filtering

```
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz
```

After filtering, kept 24 out of 24 Individuals

After filtering, kept 7486938 out of a possible 7486938 Sites

***Biallelic vs. multi-allelic SNPs*:** So, SNPs with >2 alleles probably reflect sequence or mapping errors. We also want to get rid of SNPs showing <2 alleles. ***Minor allele frequency (MAF):*** Sequencing errors are relatively common, but they tend to happen randomly and affect only 1 read at a time. Thus, if we have a SNP that is only seen very rarely, it may be a sequencing error, and should be discarded. For us, the most liberal MAF filters would be 1 allele copy out of the total 2N copies, or 1/48 = 0.02. *Missing data across individuals:* ***Missing data across individuals:*** Get rid of sites where  fewer than 80% of our samples have data. Missing data is a problem for any analysis, and population genetic statistics can behave oddly (i.e.. become biased) when a lot of individuals are missing data for a given SNP. 

```
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.8
```

After filtering, kept 5317 out of a possible 7486938 Sites

