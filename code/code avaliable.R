####### prepare ######
# lead snp           #
# gwas summary       #
# regulatory network #
######################



# Lead SNP extend
locus <- read.table("/Users/zhangjing/Desktop/hyperuricemia-main/input/locus_example.txt",header = T)
library(data.table)
trait <- fread("/Users/zhangjing/Desktop/hyperuricemia-main/input/gwassummary_example.txt",header=T)
names(trait)
trait$chr_hg38 <- paste("chr",trait$Chr,sep = '')
trait$start_hg38 <- trait$pos_hg38 -5
trait$end_hg38 <- trait$pos_hg38 +5

trait_locus <-merge(locus,trait,by="RSID")
trait_locus$ref <- "Tin"
names(trait_locus)
trait_locus <- trait_locus[,c(1,16,10,11,13:15)]
length(unique(trait_locus$RSID))

locus <- trait_locus
locus$window_left <- locus$start_hg38 - 249995
locus$window_right <- locus$end_hg38 + 249995
names(locus)
names(locus) <- c("RSID_locus","ref","Pvalue_locus","log10P_locus","chr_hg38_locus","start_hg38_locus",
                  "end_hg38_locus","window_left","window_right")
trait_locus <- locus
length(unique(trait_locus$RSID_locus))

snp <- trait
library(GenomicRanges)
snp.g=GRanges(data.frame(chrom=snp$chr_hg38,start=snp$start_hg38,end=snp$end_hg38))
trait_locus.g=GRanges(data.frame(chrom=trait_locus$chr_hg38_locus,start=trait_locus$window_left,end=trait_locus$window_right))
snp_trait_locus <- suppressWarnings(findOverlaps(trait_locus.g,snp.g))
snp_trait_locus_overlap <- data.frame(trait_locus[queryHits(snp_trait_locus),],snp[subjectHits(snp_trait_locus),])
snp_trait_locus_overlap <-snp_trait_locus_overlap[!duplicated(snp_trait_locus_overlap),]
snp_trait_locus_overlap$P.value <- as.numeric(snp_trait_locus_overlap$P.value)
snp_trait_locus_overlap <- snp_trait_locus_overlap[snp_trait_locus_overlap$P.value < 0.00000005,]
length(unique(snp_trait_locus_overlap$RSID))
length(unique(snp_trait_locus_overlap$RSID_locus))

trait_info <- snp_trait_locus_overlap
length(unique(trait_info$RSID_locus))
write.csv(trait_info,"locus_example_info.csv",row.names = F)






# LSS quantile
a <- unique(trait_info$RSID_locus)

i=1
info <- trait_info[trait_info$RSID_locus == a[i],]
cutoff <- info$log10P_locus[1]
cutoff_region <- cutoff * 0.75
info_filter <- info[info$log10P > cutoff_region,]
set <- info_filter

for (i in 2:length(a)) {
  info <- trait_info[trait_info$RSID_locus == a[i],]
  cutoff <- info$log10P_locus[1]
  cutoff_region <- cutoff * 0.75
  info_filter <- info[info$log10P > cutoff_region,]
  set <- rbind(set,info_filter)
}
length(unique(set$RSID_locus))
length(unique(locus$RSID_locus))
write.csv(set,"./locus_example_quantile.csv",row.names = F)



# Functional SNP to target gene
rm(list=ls())
trait <- read.csv("locus_example_quantile.csv",header = T)
p2g <- readRDS("kidney_p2g.rds")

snp <- trait
library(GenomicRanges)
snp.g=GRanges(data.frame(chrom=snp$chr_hg38,start=snp$start_hg38,end=snp$end_hg38))
p2g.g=GRanges(data.frame(chrom=p2g$chr,start=p2g$start,end=p2g$end))
snp_p2g <- suppressWarnings(findOverlaps(p2g.g,snp.g))
snp_p2g_overlap <- data.frame(p2g[queryHits(snp_p2g),],snp[subjectHits(snp_p2g),])
snp_p2g_overlap <-snp_p2g_overlap[!duplicated(snp_p2g_overlap),]
info_p2g <- snp_p2g_overlap
length(unique(info_p2g$RSID_locus))
length(unique(info_p2g$RSID))
length(unique(info_p2g$gene))
write.csv(info_p2g,"./locus_example_p2g.csv",row.names = F)




#GRPS
locus <- read.csv("locus_example_p2g.csv",header = T)
locus$peak <- paste(locus$chr,locus$start,locus$end,sep = '_')
names(locus)
#gene; linkage score; peak; -log10(pvalue)
locus <- locus[,c(1,4,35,30)]
locus <-locus[!duplicated(locus),]
locus$Correlation <- (locus$Correlation)^2

df <- data.frame(0,0)
names(df) <- c("gene","GRPS")
set <- df[-1,]

a <- as.character(unique(locus$gene))
info <- locus[1,]
info <- info[!duplicated(info),]
result <- info[-1,]
#highest disease association of high-risk variant in peak
for (i in 1:length(a)) {
  test <- locus[locus$gene == a[i],]
  head(test)
  peak <- unique(test$peak)
  for (j in 1:length(peak)) {
    ttest <- test[test$peak == peak[j],]
    info <- ttest[which(ttest$log10P == max(ttest$log10P)),]
    result <- rbind(result,info)
  }
  
}

length(unique(result$gene))
length(unique(locus$gene))


#calculate GRPS
locus <- result
for (i in 1:length(a)) {
  test <- locus[locus$gene == a[i],]
  head(test)
  forGRPS <- test[,c(1,2,4)]
  forGRPS$CP <- abs(forGRPS$Correlation) * forGRPS$log10P
  forGRPS <- forGRPS[,c(1,4)]
  forGRPS$CP <- round(forGRPS$CP,2)
  CPagg <- aggregate(forGRPS, list(forGRPS$gene), unique)
  
  if (dim(test)[1] == 1) {
    df$gene <- test$gene
    df$GRPS <- forGRPS$CP

  } else { 
    df$gene <- unique(test$gene)
    df$GRPS <- sum(forGRPS$CP)
  }
  
  set <- rbind(set,df)
}


length(unique(set$gene))
names(set)

set <- set[order(set$GRPS,decreasing = T),]
write.csv(set,"./locus_example_GRPS.csv",row.names = F)
