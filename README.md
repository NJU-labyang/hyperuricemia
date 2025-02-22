# hyperuricemia
A novel approach to comprehensively explore genetic contribution and nominate reliable causal genes for complex diseases

**About**

We introduce a locus-specific stratification (LSS) and gene regulatory prioritization score (GRPS) approach to comprehensively explore genetic contribution and nominate reliable causal genes for complex diseases.

**Getting Started**

This program requires R 4.0, several R packages. Versions the R packages has been tested on data.table (version = 1.16.4), GenomicRanges (version = 1.58.0)

**Demo**

We have provided example input data in the 'input' folder. You can run the example file using the code provided in the 'code' folder. The outputs generated from running above are provided in the 'output' folder.

***Input***\
gwassummary_example.txt (GWAS summary statistics)\
kidney_p2g.rds (Regulatory network)\
locus_example.txt (Lead SNP for the locus)

***output***\
locus_example_info.csv (Contains all SNPs in the locus)\
locus_example_quantile.csv (Locus-specific stratification information)\
locus_example_p2g.csv (Identification of functional SNPs)\
locus_example_GRPS.csv （GRPS information）

Expected runtime for the demo on a 'normal' desktop computer is approximately 5 to 10 minutes

**Contact**

Jing Zhang （DG20350040@smail.nju.edu.cn）
