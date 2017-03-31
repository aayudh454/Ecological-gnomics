### Page 19: 2017-03-31. Homework3_population genetics

#### 1. Filtering

```
[aadas@pbio381 homework3]$ vcftools --gzvcf SSW_by24inds.txt.vcf.gz
```

After filtering, kept 24 out of 24 Individuals

After filtering, kept 7486938 out of a possible 7486938 Sites

So, SNPs with >2 alleles probably reflect sequence or mapping errors. We also want to get rid of SNPs showing <2 alleles.