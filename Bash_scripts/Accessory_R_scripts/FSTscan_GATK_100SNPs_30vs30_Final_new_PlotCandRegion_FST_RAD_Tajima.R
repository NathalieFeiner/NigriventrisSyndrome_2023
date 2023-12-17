### This code creates Fig. 2B and D

dev.off()
rm(list=ls())
library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggthemes)
library(qqman)
library(dplyr)
library(stringr)

### FST outlier scans on WGS
Data <- read.table("FSTscan_GATK_10vs10_SA.weir.fst.averagesliding_100SNPs", header = F) #do this for both IT and SA
colnames(Data) <- c("CHROM","BP","end","FST","nSNPs")
Data$Loc <- Data$CHROM
Data$CHROM <- gsub("CM014743.1\\b", "1", Data$CHROM)
Data$CHROM <- gsub("CM014744.1\\b", "2", Data$CHROM)
Data$CHROM <- gsub("CM014745.1\\b", "3", Data$CHROM)
Data$CHROM <- gsub("CM014746.1\\b", "4", Data$CHROM)
Data$CHROM <- gsub("CM014747.1\\b", "5", Data$CHROM)
Data$CHROM <- gsub("CM014748.1\\b", "6", Data$CHROM)
Data$CHROM <- gsub("CM014749.1\\b", "7", Data$CHROM)
Data$CHROM <- gsub("CM014750.1\\b", "8", Data$CHROM)
Data$CHROM <- gsub("CM014751.1\\b", "9", Data$CHROM)
Data$CHROM <- gsub("CM014752.1\\b", "10", Data$CHROM)
Data$CHROM <- gsub("CM014753.1\\b", "11", Data$CHROM)
Data$CHROM <- gsub("CM014754.1\\b", "12", Data$CHROM)
Data$CHROM <- gsub("CM014755.1\\b", "13", Data$CHROM)
Data$CHROM <- gsub("CM014756.1\\b", "14", Data$CHROM)
Data$CHROM <- gsub("CM014757.1\\b", "15", Data$CHROM)
Data$CHROM <- gsub("CM014758.1\\b", "16", Data$CHROM)
Data$CHROM <- gsub("CM014759.1\\b", "17", Data$CHROM)
Data$CHROM <- gsub("CM014760.1\\b", "18", Data$CHROM)
Data$CHROM <- gsub("CM014761.1\\b", "19", Data$CHROM)
Data$CHROM <- as.numeric(Data$CHROM)
Data$SNP <- paste(Data$Loc,Data$BP,sep="_")
Data$zFst <- Data$FST/sd(Data$FST)
jpeg("FSTscans_GATK_100SNPs_10vs10_SA_100SNPs.jpeg", width=2500, height=300)
#pdf("FSTscans_GATK_100SNPs_10vs10_SA_100SNPs.pdf", width = 35, height = 5, useDingbats=FALSE)
manhattan(Data, chr="CHROM", p="zFst", logp = FALSE, chrlabs = c(1:18, "Z"), main ="Combined (30 vs 30)", ylab = "zFST",
          genomewideline = FALSE, suggestiveline = quantile(Data$zFst, 0.995), col = c("grey70", "grey30"))
dev.off()

#only chr. 12
Data_Chr12 <- subset(Data, CHROM==12)
pdf("../../Rscripts/FSTscans_GATK_100SNPs_10vs10_SA_100SNPs_cand.pdf", width = 15, height = 6, useDingbats=FALSE)
manhattan(Data_Chr12, chr="CHROM", p="zFst", logp = FALSE, ylab = "zFST",xlim = c(21200000, 21850000),
          genomewideline = FALSE, suggestiveline = quantile(Data$zFst, 0.995), col = c("grey30"))
dev.off()

#jpeg("../../Rscripts/FSTscans_GATK_new_SA_100SNPs.jpeg", width=3000, height=1600)
#pdf("../../Rscripts/FSTscans_GATK_100SNPs_30vs30_Chr12_new_100kb25.pdf", width = 14  , height = 8, useDingbats=FALSE)
  p1 <- ggplot(Data_Chr12, aes(x=BP, y=zFst)) +
  geom_point(size=2) + xlim(20300000,22800000) +
  scale_color_manual(values = "darkgrey") +
  theme_bw() + theme(legend.position="none", plot.title = element_text(size=20)) + ggtitle("FST-scan") +
  geom_hline(yintercept=quantile(Data$zFst, 0.995), linetype="dashed", color = "red")
#dev.off()

### GWAS on RAD-seq data
Plink_adjusted <- read.table("/IT_PC/plink.assoc.linear.adjusted", header=T) #do this for both IT and SA
Plink_nonadjusted <- read.table("/IT_PC/plink.assoc.linear", header=T)
Plink_nonadjusted <- Plink_nonadjusted[,1:3]
Plink_nonadjusted <- Plink_nonadjusted[!duplicated(Plink_nonadjusted), ]
Plink_new <- merge(Plink_adjusted, Plink_nonadjusted, by.x = "SNP", by.y = "SNP", all.x=T, sort=F)
Plink_new <- Plink_new[,c(2,1,12,3:10)]
Plink_new$logP <- -log(Plink_new$FDR_BH)
Plink_new$Loc <- Plink_new$CHR.x
Plink_new$CHR.x <- gsub("CM014743.1\\b", "1", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014744.1\\b", "2", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014745.1\\b", "3", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014746.1\\b", "4", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014747.1\\b", "5", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014748.1\\b", "6", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014749.1\\b", "7", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014750.1\\b", "8", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014751.1\\b", "9", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014752.1\\b", "10", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014753.1\\b", "11", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014754.1\\b", "12", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014755.1\\b", "13", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014756.1\\b", "14", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014757.1\\b", "15", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014758.1\\b", "16", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014759.1\\b", "17", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014760.1\\b", "18", Plink_new$CHR.x)
Plink_new$CHR.x <- gsub("CM014761.1\\b", "19", Plink_new$CHR.x)
Plink_new$CHR.x <- as.numeric(Plink_new$CHR.x)
jpeg("FSTscans_RADseq_IT.jpeg", width=2500, height=300)
#pdf("FSTscans_RADseq_IT.pdf", width = 35, height = 5, useDingbats=FALSE)
manhattan(Plink_new, chr="CHR.x", p="FDR_BH", chrlabs = c(1:18, "Z"), main ="RADseq_SA",
          genomewideline = FALSE, suggestiveline = -log10(0.05), col = c("grey70", "grey30"))
dev.off()

### delta Tajima's D on WGS data
Data_brown <- read.table("/crex/proj/snic2020-6-73/Projects/Wallie_GWAS/Analysis_149samples/GATK/TajimasD_pi/Final_new/WG_5kb_TajimaD_brown_IT.Tajima.D", header = T) #do this for both IT and SA
Data_green <- read.table("/crex/proj/snic2020-6-73/Projects/Wallie_GWAS/Analysis_149samples/GATK/TajimasD_pi/Final_new/WG_5kb_TajimaD_green_IT.Tajima.D", header = T) #do this for both IT and SA
Data_brown$POS <- paste(Data_brown$CHROM,Data_brown$BIN_START, sep="_")
Data_green$POS <- paste(Data_green$CHROM,Data_green$BIN_START, sep="_")
Data <- merge(Data_green, Data_brown, by.x="POS", by.y="POS")
Data <- Data[,-c(2,3,4,6:8)]
colnames(Data) <- c("POS", "green", "brown")
Data$DeltaTajima <- Data$brown-Data$green
Data[c('CHROM', 'BP')] <- str_split_fixed(Data$POS, '_', 2)
Data$Loc <- Data$CHROM
Data$CHROM <- gsub("CM014743.1\\b", "1", Data$CHROM)
Data$CHROM <- gsub("CM014744.1\\b", "2", Data$CHROM)
Data$CHROM <- gsub("CM014745.1\\b", "3", Data$CHROM)
Data$CHROM <- gsub("CM014746.1\\b", "4", Data$CHROM)
Data$CHROM <- gsub("CM014747.1\\b", "5", Data$CHROM)
Data$CHROM <- gsub("CM014748.1\\b", "6", Data$CHROM)
Data$CHROM <- gsub("CM014749.1\\b", "7", Data$CHROM)
Data$CHROM <- gsub("CM014750.1\\b", "8", Data$CHROM)
Data$CHROM <- gsub("CM014751.1\\b", "9", Data$CHROM)
Data$CHROM <- gsub("CM014752.1\\b", "10", Data$CHROM)
Data$CHROM <- gsub("CM014753.1\\b", "11", Data$CHROM)
Data$CHROM <- gsub("CM014754.1\\b", "12", Data$CHROM)
Data$CHROM <- gsub("CM014755.1\\b", "13", Data$CHROM)
Data$CHROM <- gsub("CM014756.1\\b", "14", Data$CHROM)
Data$CHROM <- gsub("CM014757.1\\b", "15", Data$CHROM)
Data$CHROM <- gsub("CM014758.1\\b", "16", Data$CHROM)
Data$CHROM <- gsub("CM014759.1\\b", "17", Data$CHROM)
Data$CHROM <- gsub("CM014760.1\\b", "18", Data$CHROM)
Data$CHROM <- gsub("CM014761.1\\b", "19", Data$CHROM)
Data$CHROM <- as.numeric(Data$CHROM)
Data$BP <- as.numeric(Data$BP)
Data$SNP <- paste(Data$Loc,Data$BP,sep="_")
Data <- Data[complete.cases(Data[ , 4]),]
jpeg("FSTscans_Tajima_IT.jpeg", width=2500, height=300)
#pdf("FSTscans_Tajima_SA.pdf", width = 35, height = 5, useDingbats=FALSE)
manhattan(Data, chr="CHROM", p="DeltaTajima", logp = F, chrlabs = c(1:18, "Z"), main ="TajimasD IT",
          genomewideline = FALSE, suggestiveline = quantile(Data$DeltaTajima, 0.995, na.rm = T), col = c("grey70", "grey30"))
dev.off()
