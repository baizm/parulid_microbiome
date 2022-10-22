setwd("~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision")
#library(ctc)
library(vegan)
library(scales)
library(ape)
library(phytools)
load('mantel_cophenetic2.RData')

ht<-read.tree('~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/eliot_example/scaledWarbler_full.tre')
plot(ht)

#convert tip labels to banding codes
ht$tip.label<-c('OVEN','WEWA','BAWW','NAWA','COYE','HOWA','AMRE','BTBW','NOPA','MAWA','BLBW','CSWA','MYWA','BTNW','CAWA')
coph<-cophenetic(ht)
coph<-as.dist(coph, diag=T, upper=F) #dist object for hclust

bray_full<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/bray_full.rds')

#use coph values
#match sample id with species names for cat distances
m<-read.csv('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/metadata_16s_merged.csv', header=T)
m<-m[order(m$id),]
rownames(m)<-m$id
m<-m[,2:31] #get rid of id column
m<-m[9:412,] #and the negatives
m2<-m[which(rownames(m) %in% attr(bray_full, 'Labels')),] #and samples not included

ids<-data.frame(id=rownames(m2), species=m2$Species)
#make sure in same order
identical(attr(bray_full, 'Labels'), ids$id) #true

l<-list('rownames'=ids$species, 'colnames'=ids$species)
mat<-matrix(nrow=270, ncol=270, dimnames=l) #empty matrix with species names as row/col names in order of m$id, should be same order as distance matrix

coph2<-as.matrix(coph)

#populate matrix with values from cat distances
mat2<-coph2[rownames(mat), colnames(mat)]
mat3<-as.dist(mat2, diag=F, upper=F) #dist object for mantel


#import distances matrices from phyloseq session (individual level matrix)
#####-----gm-cophenetic distance--------
###full dataset
#bray-curtis
identical(attr(bray_full, "Labels"), attr(mat3, "Labels")) #but its okay bc one is species codes and one is ids, but they match!!
mantel(mat3, bray_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.02, p=0.21

jac_full<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/jac_full.rds')
mantel(mat3, jac_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.04, p=0.067

uni_full<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/unifrac_full.rds')
mantel(mat3, uni_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.10, p=0.005

wuni_full<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/wunifrac_full.rds')
mantel(mat3, wuni_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.01, p=0.421

###batch 2, diet-gm
bray_b2<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/bray_b2.rds')

m_b2<-read.csv('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/metadata_16s_merged.csv', header=T)
m_b2<-m_b2[order(m_b2$id),]
rownames(m_b2)<-m_b2$id
m_b2<-m_b2[,2:31] #get rid of id column
m_b2<-m_b2[9:412,] #and the negatives
m2_b2<-m_b2[which(rownames(m_b2) %in% attr(bray_b2, 'Labels')),] #and samples not included

ids_b2<-data.frame(id=rownames(m2_b2), species=m2_b2$Species)
#make sure in same order
identical(attr(bray_b2, 'Labels'), ids_b2$id) #true

l_b2<-list('rownames'=ids_b2$species, 'colnames'=ids_b2$species)
mat_b2<-matrix(nrow=111, ncol=111, dimnames=l_b2) #empty matrix with species names as row/col names in order of m$id, should be same order as distance matrix

#coph2<-as.matrix(coph)

#populate matrix with values from cat distances
mat2_b2<-coph2[rownames(mat_b2), colnames(mat_b2)]
mat3_b2<-as.dist(mat2_b2, diag=F, upper=F) #dist object for mantel

identical(attr(bray_b2, "Labels"), attr(mat3_b2, "Labels")) #F, but check  head(ids_b2)--good!
mantel(mat3_b2, bray_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.07, p=0.112

jac_b2<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/jac_b2.rds')
mantel(mat3_b2, jac_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.03, p=0.328

unifrac_b2<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/unifrac_b2.rds')
mantel(mat3_b2, unifrac_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.04, p=0.236

wunifrac_b2<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_merged/wunifrac_b2.rds')
mantel(mat3_b2, wunifrac_b2, method = "spearman", permutations = 999, na.rm = TRUE) #-0.02, p=0.643

###batch 1
bray_b1<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/bray_b1.rds')

m_b1<-read.csv('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/metadata2_16s_ny.tsv', header=T, sep='\t')
m_b1<-m_b1[order(m_b1$id),]
rownames(m_b1)<-m_b1$id
m_b1<-m_b1[,2:29] #get rid of id column
m2_b1<-m_b1[which(rownames(m_b1) %in% attr(bray_b1, 'Labels')),] #and samples not included ##161 is correct number
identical(attr(bray_b1, 'Labels'), rownames(m2_b1)) #yes

ids_b1<-data.frame(id=rownames(m2_b1), species=m2_b1$Species)

l_b1<-list('rownames'=ids_b1$species, 'colnames'=ids_b1$species)
mat_b1<-matrix(nrow=161, ncol=161, dimnames=l_b1) #empty matrix with species names as row/col names in order of m$id, should be same order as distance matrix

#populate matrix with values from cat distances
#get rid of worm-eating, not in batch 1
coph2_b1<-coph2[-2,]
coph2_b1<-coph2_b1[,-2]

#populate matrix with values from coph distances
mat2_b1<-coph2_b1[rownames(mat_b1), colnames(mat_b1)]
mat3_b1<-as.dist(mat2_b1, diag=F, upper=F) #dist object for mantel

identical(attr(bray_b1, "Labels"), attr(mat3_b1, "Labels")) #F, but check  head(ids_b1)--good!
mantel(mat3_b1, bray_b1, method = "spearman", permutations = 999, na.rm = TRUE) #0.09, p=0.025

jac_b1<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/jac_b1.rds')
mantel(mat3_b1, jac_b1, method = "spearman", permutations = 999, na.rm = TRUE) #0.18, p=0.001

unifrac_b1<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/unifrac_b1.rds')
mantel(mat3_b1, unifrac_b1, method = "spearman", permutations = 999, na.rm = TRUE) #0.19, p=0.002

wunifrac_b1<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/wunifrac_b1.rds')
mantel(mat3_b1, wunifrac_b1, method = "spearman", permutations = 999, na.rm = TRUE) #-0.0003, p=0.506

###batch 1--NY
bray_ny<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/bray_ny.rds')

m_ny<-read.csv('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/metadata2_16s_ny.tsv', header=T, sep='\t')
m_ny<-m_ny[order(m_ny$id),]
rownames(m_ny)<-m_ny$id
m_ny<-m_ny[,2:29] #get rid of id column
m2_ny<-m_ny[which(rownames(m_ny) %in% attr(bray_ny, 'Labels')),] #and samples not included ##127 is correct number
identical(attr(bray_ny, 'Labels'), rownames(m2_ny)) #yes

ids_ny<-data.frame(id=rownames(m2_ny), species=m2_ny$Species)

l_ny<-list('rownames'=ids_ny$species, 'colnames'=ids_ny$species)
mat_ny<-matrix(nrow=127, ncol=127, dimnames=l_ny) #empty matrix with species names as row/col names in order of m$id, should be same order as distance matrix

#populate matrix with values from coph distances
mat2_ny<-coph2_b1[rownames(mat_ny), colnames(mat_ny)] #use coph2_b1, WEWA removed already
mat3_ny<-as.dist(mat2_ny, diag=F, upper=F) #dist object for mantel

identical(attr(bray_ny, "Labels"), attr(mat3_ny, "Labels")) #F, but check  head(ids_ny)--good!
mantel(mat3_ny, bray_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.15, p=0.004

jac_ny<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/jac_ny.rds')
mantel(mat3_ny, jac_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.23, p=0.001

unifrac_ny<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/unifrac_ny.rds')
mantel(mat3_ny, unifrac_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.27, p=0.001

wunifrac_ny<-readRDS('/Users/marcellabaiz/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/wunifrac_ny.rds')
mantel(mat3_ny, wunifrac_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.07, p=0.15

###------diet-cophenetic-------
#in batch 1
bray_diet2<-readRDS('~/Documents/REU2021/warbler_data/bray_diet2.rds')

m2_diet_b1<-m_b1[which(rownames(m_b1) %in% attr(bray_diet2, 'Labels')),] #and samples not included ##130 is correct number
identical(attr(bray_diet2, 'Labels'), rownames(m2_diet_b1)) #yes

ids_diet_b1<-data.frame(id=rownames(m2_diet_b1), species=m2_diet_b1$Species)

l_diet_b1<-list('rownames'=ids_diet_b1$species, 'colnames'=ids_diet_b1$species)
mat_diet_b1<-matrix(nrow=130, ncol=130, dimnames=l_diet_b1) #empty matrix with species names as row/col names in order of m$id, should be same order as distance matrix

#populate matrix with values from coph distances
mat2_diet_b1<-coph2_b1[rownames(mat_diet_b1), colnames(mat_diet_b1)] #use coph2_b1, no WEWA
mat3_diet_b1<-as.dist(mat2_diet_b1, diag=F, upper=F) #dist object for mantel

identical(attr(bray_diet2, "Labels"), attr(mat3_diet_b1, "Labels")) #F, but check  head(ids_diet_b1)--good!
mantel(mat3_diet_b1, bray_diet2, method = "spearman", permutations = 999, na.rm = TRUE) #0.08, p=0.073

jac_diet<-readRDS('~/Documents/REU2021/warbler_data/jac_diet.rds')
mantel(mat3_diet_b1, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.1014, p=0.027

unifrac_diet<-readRDS('~/Documents/REU2021/warbler_data/unifrac_diet.rds')
mantel(mat3_diet_b1, unifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.13, p=0.004

weightedu_diet<-readRDS('~/Documents/REU2021/warbler_data/weightedu_diet.rds')
mantel(mat3_diet_b1, weightedu_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.04, p=0.196

####---plots------
#batch 1 gm-cophenetic
pdf(file='figs/mantel_bc_gm-phylo3.pdf', height=4, width=4.5)
plot(as.vector(bray_b1), as.vector(mat3_b1), pch = 16, col = alpha('sandybrown', 0.4),          
     xlab='Microbiome distance (Bray-Curtis)',
     ylab='Evolutionary distance')
abline(lm(as.vector(mat3_b1)~as.vector(bray_b1)), xpd=F, lty=2, lwd=1.5, col='black')
dev.off()
#batch 1 diet-cophenetic
pdf(file='figs/mantel_bc_diet-phylo3.pdf', height=4, width=4.5)
plot(as.vector(bray_diet2),as.vector(mat3_diet_b1), pch = 16, col = alpha('aquamarine4',0.4),          
     xlab='Diet distance (Bray-Curtis)',
     ylab='Evolutionary distance')
abline(lm(as.vector(mat3_diet_b1)~as.vector(bray_diet2)), xpd=F, lty=2, lwd=1.5, col='black')
dev.off()

###-----iet-gm mantels--------
#diet-gm, full dataset
bray_micro_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/bray_micro_full.rds')
bray_diet_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/bray_diet_full.rds')
identical(attr(bray_micro_full, "Labels"), attr(bray_diet_full, "Labels")) #yes
mantel(bray_diet_full, bray_micro_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.05884, p=0.046

jac_micro_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/jac_micro_full.rds')
jac_diet_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/jac_diet_full.rds')
mantel(jac_diet_full, jac_micro_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.1744, p=0.001

unifrac_micro_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/unifrac_micro_full.rds')
unifrac_diet_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/unifrac_diet_full.rds')
mantel(unifrac_diet_full, unifrac_micro_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.1632, p=0.001

wunifrac_micro_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/wunifrac_micro_full.rds')
wunifrac_diet_full<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/wunifrac_diet_full.rds')
mantel(wunifrac_diet_full, wunifrac_micro_full, method = "spearman", permutations = 999, na.rm = TRUE) #0.03, p=0.324

#diet-gm batch 2
bray_micro_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/bray_micro_b2.rds')
bray_diet_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/bray_diet_b2.rds')
mantel(bray_diet_b2, bray_micro_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.008, p=0.426

jac_micro_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/jac_micro_b2.rds')
jac_diet_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/jac_diet_b2.rds')
mantel(jac_diet_b2, jac_micro_b2, method = "spearman", permutations = 999, na.rm = TRUE) #-0.03, p=0.684

unifrac_micro_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/unifrac_micro_b2.rds')
unifrac_diet_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/unifrac_diet_b2.rds')
mantel(unifrac_diet_b2, unifrac_micro_b2, method = "spearman", permutations = 999, na.rm = TRUE) #-0.016, p=0.577

wunifrac_micro_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/wunifrac_micro_b2.rds')
wunifrac_diet_b2<-readRDS('~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/wunifrac_diet_b2.rds')
mantel(wunifrac_diet_b2, wunifrac_micro_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.056, p=0.247

#diet-gm batch 1
bray_micro2<-readRDS("~/Documents/REU2021/warbler_data/bray_micro2.rds")
bray_diet2<-readRDS("~/Documents/REU2021/warbler_data/bray_diet2.rds")
mantel(bray_micro2, bray_diet2, method = "spearman", permutations = 999, na.rm = TRUE) #0.09, p=0.046

jac_micro<-readRDS("~/Documents/REU2021/warbler_data/jac_micro.rds")
jac_diet<-readRDS("~/Documents/REU2021/warbler_data/jac_diet.rds")
mantel(jac_micro, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.213, p=0.001

unifrac_micro<-readRDS("~/Documents/REU2021/warbler_data/unifrac_micro.rds")
unifrac_diet<-readRDS("~/Documents/REU2021/warbler_data/unifrac_diet.rds")
mantel(unifrac_micro, unifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.22, p=0.001

weightedu_micro<-readRDS("~/Documents/REU2021/warbler_data/weightedu_micro.rds")
weightedu_diet<-readRDS("~/Documents/REU2021/warbler_data/weightedu_diet.rds")
mantel(weightedu_micro, weightedu_diet, method = "spearman", permutations = 999, na.rm = TRUE) #-0.04771, p=0.774

#diet-gm, batch 1-NY
bray_micro_ny<-readRDS("~/Documents/REU2021/warbler_data/bray_micro_ny.rds")
bray_diet_ny<-readRDS("~/Documents/REU2021/warbler_data/bray_diet_ny.rds")
mantel(bray_micro_ny, bray_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) 

jac_micro_ny<-readRDS("~/Documents/REU2021/warbler_data/jac_micro_ny.rds")
jac_diet_ny<-readRDS("~/Documents/REU2021/warbler_data/jac_diet_ny.rds")
mantel(jac_micro_ny, jac_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) 

uni_micro_ny<-readRDS("~/Documents/REU2021/warbler_data/uni_micro_ny.rds")
uni_diet_ny<-readRDS("~/Documents/REU2021/warbler_data/uni_diet_ny.rds")
mantel(uni_micro_ny, uni_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) 

wuni_micro_ny<-readRDS("~/Documents/REU2021/warbler_data/wuni_micro_ny")
wuni_diet_ny<-readRDS("~/Documents/REU2021/warbler_data/wuni_diet_ny.rds")
mantel(wuni_micro_ny, wuni_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) 

###----plot----
#batch 1 gm-diet
pdf(file='figs/mantel_bc_gm-diet2.pdf', height=4, width=4.5)
plot(as.vector(bray_micro2), as.vector(bray_diet2), pch = 16, col = alpha('gray50',0.4),          
     xlab='Microbiome distance (Bray-Curtis)',
     ylab='Diet distance (Bray-Curtis)')
abline(lm(as.vector(bray_diet2)~as.vector(bray_micro2)), xpd=F, lty=2, lwd=1.5, col='black')
dev.off()

save.image('mantel_cophenetic2.RData')

#check out other plots for fun
plot(as.vector(jac_micro), as.vector(jac_diet), pch = 16, col = alpha('gray50',0.4),          
     xlab='Microbiome distance (Jaccard)',
     ylab='Diet distance (Jaccard)', xlim=c(0.55,1), ylim=c(0.55,1))
abline(lm(as.vector(jac_diet)~as.vector(jac_micro)), xpd=F, lty=2, lwd=1.5, col='black')

plot(as.vector(unifrac_micro), as.vector(unifrac_diet), pch = 16, col = alpha('gray50',0.4),          
     xlab='Microbiome distance (Unweighted UniFrac)',
     ylab='Diet distance (Unweighted UniFrac)')
     #xlim=c(0.55,1), ylim=c(0.55,1))
abline(lm(as.vector(unifrac_diet)~as.vector(unifrac_micro)), xpd=F, lty=2, lwd=1.5, col='black')

plot(as.vector(weightedu_micro), as.vector(weightedu_diet), pch = 16, col = alpha('gray50',0.4),          
     xlab='Microbiome distance (Weighted UniFrac)',
     ylab='Diet distance (Weighted UniFrac)')
#, xlim=c(0.55,1), ylim=c(0.55,1))
abline(lm(as.vector(weightedu_diet)~as.vector(weightedu_micro)), xpd=F, lty=2, lwd=1.5, col='black')
