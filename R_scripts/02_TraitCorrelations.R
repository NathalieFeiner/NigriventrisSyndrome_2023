rm(list=ls())

library(openxlsx)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(pcaMethods)
library(dplyr)

options(scipen=999)

setwd("C:/Users/Nathalie/Dropbox/WallLizard_NanoporeAssemblies/GreenPackage/PhenotypeDataAnalyses")

#Import the current new data set - this includes data up till 2018
data <- read.csv("Phenotype_data_final_June2023.csv")

#Import lineage assignments per population (from structure on RAD-seq data plus extrapolation)
Pops <- read.csv("PopsLocs.csv")
Pops <- Pops[,-c(3,4)]
data <- merge(data, Pops, by = "abbpop", all.x=T)

#Greenness per pop
Green_mean_sex <- as.data.frame(data %>% group_by(abbpop, sex) %>% dplyr::summarize(Green = mean(GreenResolved, na.rm=TRUE)))
Green_mean_sex <- subset(Green_mean_sex, sex=="M")
Green_mean_sex <- Green_mean_sex[,-c(2)]
#Green_mean <- data %>% group_by(abbpop) %>% dplyr::summarize(Green = mean(GreenResolved, na.rm=TRUE))
data <- merge(data, Green_mean_sex, by="abbpop")

Pops <- merge(Pops, Green_mean_sex, by="abbpop")
Pops_SA_intro <- subset(Pops, lineage=="SA" & Green>=3)
Pops_SA_intro$lineage_intro <- "SA_intro"
Pops_SA_intro <- Pops_SA_intro[,-c(2,3)]

data <- merge(data,Pops_SA_intro, by="abbpop", all.x=T)

#change some colnames
colnames(data)[c(11,14,15,16)] <- c("green","black","blue","rel.headlength")

#subset data by lineage
data_IT_M <- subset(data, lineage=="ITA" & sex=="M")
data_Hyb_M <- subset(data, lineage=="MIXED" & sex=="M")
data_SAintro_M <- subset(data, lineage_intro == "SA_intro" & sex=="M")
data_IT_F <- subset(data, lineage=="ITA" & sex=="F")
data_Hyb_F <- subset(data, lineage=="MIXED" & sex=="F")
data_SAintro_F <- subset(data, lineage_intro == "SA_intro" & sex=="F")

#correlation IT males
matrix_IT_M <- as.matrix(scale(data_IT_M[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_M, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("IT_M_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation hybrid males
matrix_IT_M <- as.matrix(scale(data_Hyb_M[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_M, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("Hyb_M_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation SAintro males
matrix_IT_M <- as.matrix(scale(data_SAintro_M[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_M, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("SAintro_M_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation IT females
matrix_IT_F <- as.matrix(scale(data_IT_F[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_F, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("IT_F_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation hybrid females
matrix_Mix_F <- as.matrix(scale(data_Hyb_F[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_Mix_F, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("Hyb_F_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation SAintro females
matrix_IT_F <- as.matrix(scale(data_SAintro_F[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_F, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("SAintro_F_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

###Phenotypic integration score
library(geomorph)

#compare integration in males between IT, Hyb and SAintro
data_M <- subset(data, (lineage_intro=="SA_intro" | lineage== "ITA" | lineage== "MIXED") & sex=="M")

#delete incomplete obs
data_M <- data_M[complete.cases(data_M[,c(8,16,11,14,15)]), ]

###impute missing data
#data_sub <- data_F[,c(8,16,11,14,15)]
#summary(data_sub)
#pc <- pca(data_sub, nPcs=5, method="ppca", seed=123)
#imputed <- as.data.frame(completeObs(pc))
#data_M <- data_M[,-c(8,16,11,14,15)]
#data_F <- data_F[,-c(8,16,11,14,15)]
#data_final <- cbind(data_M, imputed)
#data_final <- cbind(data_F, imputed)

data_M[,c(8,16,11,14,15)] <- scale(data_M[,c(8,16,11,14,15)])
data_M_IT <- subset(data_M, lineage=="ITA")
data_M_Hyb <- subset(data_M, lineage=="MIXED")
data_M_SAintro <- subset(data_M, lineage_intro=="SA_intro")

IT_M <- integration.Vrel(as.matrix(data_M_IT[,c(8,16,11,14,15)])) #extract $Re.obs as Vrel
HY_M <- integration.Vrel(as.matrix(data_M_Hyb[,c(8,16,11,14,15)])) #extract $Re.obs as Vrel
SA_M <- integration.Vrel(as.matrix(data_M_SAintro[,c(8,16,11,14,15)])) #extract $Re.obs as Vrel
res <- compare.ZVrel(IT_M, HY_M, SA_M) #extract P-value for difference

#same for females
data_F <- subset(data, (lineage_intro=="SA_intro" | lineage== "ITA" | lineage== "MIXED") & sex=="F")

data_F <- data_F[complete.cases(data_F[,c(8,16,11,14,15)]), ]

data_F[,c(8,16,11,14,15)] <- scale(data_F[,c(8,16,11,14,15)])
data_F_IT <- subset(data_F, lineage=="ITA")
data_F_Hyb <- subset(data_F, lineage=="MIXED")
data_F_SAintro <- subset(data_F, lineage_intro=="SA_intro")

IT_F <- integration.Vrel(as.matrix(data_F_IT[,c(8,16,11,14,15)]))
HY_F <- integration.Vrel(as.matrix(data_F_Hyb[,c(8,16,11,14,15)]))
SA_F <- integration.Vrel(as.matrix(data_F_SAintro[,c(8,16,11,14,15)]))
res <- compare.ZVrel(IT_F, HY_F, SA_F)


_________________________________
delete everything from here on!

# calculate average for each population and for males and females separately
mphen<-aggregate(cbind(mmass=data$mass, mrhl=data$rel.headlength, mgreen=data$green, mblack=data$black, mblue=data$blue),
                 by=list(abbpop=data$abbpop, sex=data$sex), mean,na.rm=TRUE)

mphen <- merge(mphen,Pops, by = "abbpop")

mphenM_ITA <- subset(mphen, sex =="M" & lineage =="ITA")
mphenF_ITA <- subset(mphen, sex =="F" & lineage == "ITA")

MmatrixITA <- as.matrix(scale(mphenM_ITA[,c(3:7)]))
cor(MmatrixITA, method = "pearson", use = "complete.obs")

ITA_M_cor <- rcorr(MmatrixITA)
ITA_M_cor

plot(mphenM_ITA$msvl,mphenM_ITA$mmass)
plot(mphenM_ITA$msvl,mphenM_ITA$mgreen)
plot(mphenM_ITA$msvl,mphenM_ITA$mblack)
plot(mphenM_ITA$msvl,mphenM_ITA$mblue)
plot(mphenM_ITA$mblack,mphenM_ITA$mgreen)
plot(mphenM_ITA$mblack,mphenM_ITA$mblue)
plot(mphenM_ITA$mgreen,mphenM_ITA$mblue)

plot(mphenM_ITA$mblack,mphenM_ITA$mrhl)







----
  
  # Old stuff for the different clines
  
  clineApu <- subset(data2012_2018, abbpop == "FDM" | abbpop == "SZ" | abbpop == "CSA" | abbpop == "RN" | abbpop == "LV" | abbpop == "LVF" | abbpop == "VF") # Apuane cline
clineApu$mountain <- "Apuane"

clineApp <- subset(data2012_2018, abbpop == "STM" | abbpop == "CHI" | abbpop == "BLI" | abbpop == "BZ" |abbpop == "CTG" | abbpop == "PZ" | abbpop == "PSO" | abbpop == "BVN" | abbpop == "MOI" | abbpop == "MOB") #Appenine cline
clineApp$mountain <- "Appennine"

# cline along coast
clinecoast <- subset(data2012_2018, abbpop == "QC" | abbpop == "VI" | abbpop == "ST" |abbpop == "LE" | abbpop == "SL" | abbpop == "RA" | abbpop == "GN" | abbpop == "ME" | abbpop == "VA" | abbpop == "NL" | abbpop == "LO") # coast cline

cline_loc <- data.frame(abbpop = c("QC","VI","ST","LE","SL","RA","GN","ME","VA","NL","LO"), pop_order = c("A_QC","B_VI","C_ST","D_LE","E_SL","F_RA","G_GN","H_ME","I_VA","J_NL","K_LO"))

clinecoast <- merge(cline_loc,clinecoast, by ="abbpop") 

attach(clinecoast)
clinecoast <- clinecoast[order(pop_order),]
detach(clinecoast)

# Rome cline
clineRome <- subset(data2012_2018, abbpop == "FU" | abbpop == "RO" | abbpop == "LC" |abbpop == "PS" | abbpop == "MF" | abbpop == "OV" | abbpop == "PX" | abbpop == "LS" | abbpop == "FG" | abbpop == "DS" | abbpop == "AS") # Rome cline

Rome_loc <- data.frame(abbpop = c("FU","RO","LC","PS","MF", "OV","PX","LS","FG","DS","AS"), pop_order = c("A_FU","B_RO","C_LC","D_PS","E_MF", "F_OV","G_PX","H_LS","I_FG","J_DS","K_AS"))

clineRome <- merge(Rome_loc,clineRome, by ="abbpop") 

attach(clineRome)
clineRome <- clineRome[order(pop_order),]
detach(clineRome)



```{r}


ApuM <- subset(colourApu, sex =="M")
AppM <- subset(colourApp, sex =="M")

ggplot(bothclines, aes(x=svl, y=GeoffMara, color=sex)) +
  geom_point()

cor.test(AppM$svl,AppM$GreenResolved)


#old crap for analyses of phenotype trajectories

library(RRPP)

colourApp <- clineApp[,c(1,2,5,6,11,16,21,22,29,30)]

colourApp <- colourApp[complete.cases(colourApp), ]

colourApu <- clineApu[,c(1,2,5,6,11,16,21,22,29,30)]

colourApu <- colourApu[complete.cases(colourApu), ]

colourboth <- bothclines[,c(1,2,5,6,11,16,21,22,29,30)]

colourboth <- colourboth[complete.cases(colourboth), ]

colourboth_M <- subset(colourboth, sex =="M")

Appmatrix <- as.matrix(scale(colourApp[,c(5,6,9)]))

Apumatrix <- as.matrix(scale(colourApu[,c(5,6,9)]))

bothmatrix <- as.matrix(scale(colourboth[,c(5,6,9)]))

colourcoast <- clinecoast[,c(2,3,6:12,17,28)]

colourcoast <- colourcoast[complete.cases(colourcoast),]
coastmatrix <- as.matrix(scale(colourcoast[,c(9,10,11)]))

coastmatrix_bg <- as.matrix(scale(colourcoast[,c(9,10)]))

colourRome <- clineRome[,c(2,3,6:12,17,28)]

colourRome <- colourRome[complete.cases(colourRome),]
Romematrix <- as.matrix(scale(colourRome[,c(9,10,11)]))

Romematrix_bg <- as.matrix(scale(colourRome[,c(9,10)]))

#trajectory <- as.vector(c("QC","VI","MG","ST","LE","SL","RA","GN","ME","VA","NL","LO"))

#bothmatrix_M <- as.matrix(scale(colourboth_M[,c(5,6)]))

####################

fit <- lm.rrpp(Romematrix_bg ~ svl + pop_order*sex, data = colourRome, iter = 999)
reveal.model.designs(fit)
TA <- trajectory.analysis(fit, groups = colourRome$sex,
                          traj.pts = colourRome$pop_order, print.progress = FALSE)
summary(TA, attribute = "MD", show.trajectories=TRUE) # Magnitude difference (absolute difference between path distances)
summary(TA, attribute = "TC", angle.type = "deg", show.trajectories=TRUE) # Correlations (angles) between trajectories
summary(TA, attribute = "SD", show.trajectories=TRUE) # No shape differences between vectors
# Retain results
TA.summary <- summary(TA, attribute = "SD")
TA.summary$summary.table
# Plot results
TP <- plot(TA, pch = as.numeric(colourcoast$pop_order) + 20, bg = as.numeric(colourcoast$sex),
           cex = 0.7, col = "gray")
add.trajectories(TP, traj.pch = c(21, 22), start.bg = 1, end.bg = 11)
legend("topright", levels(colourcoast$abbpop), pch = c(21, 22), pt.bg = 1)

summary.trajectory.analysis(TA, attribute ="MD")

####################

fit <- lm.rrpp(Appmatrix ~ svl + sex*alt, SS.type = "III",
               data = colourApp, print.progress = FALSE, iter = 999)
summary(fit, formula = TRUE)
anova(fit)
coef(fit, test = TRUE)

# Predictions (holding alternative effects constant)
shapeDF <- data.frame(sex = c("M","M","M","M","M","M","M","M","M","M","M","F", "F", "F", "F", "F", "F", "F","F", "F", "F", "F"), pop_order = c("A_QC","B_VI","C_ST","D_LE","E_SL","F_RA","G_GN","H_ME","I_VA","J_NL","K_LO","A_QC","B_VI","C_ST","D_LE","E_SL","F_RA","G_GN","H_ME","I_VA","J_NL","K_LO")) # levels(colourApp$sex), alt = 

shapeDF <- data.frame(sex = c("M","M","M","M","M","M","M","M","M","M","M", "F", "F", "F", "F", "F", "F", "F","F", "F", "F", "F"), pop_order = c("A_FU","B_RO","C_LC","D_PS","E_MF","F_OV", "G_PX","H_LS","I_FG","J_DS","K_AS","A_FU","B_RO","C_LC","D_PS","E_MF","F_OV", "G_PX","H_LS","I_FG","J_DS","K_AS")) # 

shapeDF <- data.frame(sex = c("M","M","M","M","F","F","F","F"), alt = c("a0","a250","a500","a1000","a0","a250","a500","a1000"))
#rownames(shapeDF) <- paste(shapeDF$sex, shapeDF$alt, sep = ".")
shapeDF
shapePreds <- predict(fit, shapeDF, confidence = 0.95)
summary(shapePreds)
summary(shapePreds, PC = TRUE)
shapePreds$pca
# Plot prediction
plot(shapePreds, PC = TRUE)
plot(shapePreds, PC = TRUE, ellipse = TRUE)

plot(shapePreds, PC = TRUE, ellipse = TRUE,
     pch = 19, col = 1:NROW(shapeDF))

plot(shapePreds, PC = TRUE, ellipse = TRUE,
     pch = 19, col = as.rownames(groups))
# Diagnostics plots of residuals
plot(fit)
# PC-plot of fitted values
groups <- interaction(colourRome$sex, colourRome$pop_order)
plot(fit, type = "PC", pch = 19, col = as.numeric(groups))
# Regression-like plot
plot(fit, type = "regression", reg.type = "PredLine",
     predictor = colourApp$svl, pch=19,
     col = as.numeric(groups))