setwd("~/Documents/Toews_Lab/warbler_diet/diet_and_microbiome_merged/")
library(phyloseq)
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(reshape2)
library(decontam)
library(gridExtra)
library(plyr)
library(Rmisc)

#read in raw diet physeq object (2017-2020 data)
ps_diet_raw<-readRDS(file = "../diet_raw_physeq.rds")

#read in filtered and rarefied 16s data (2017-2020 data)
########-------create physeq object---------------------########
####------read in OTU table, has to be a matrix---------------------
otu_table<-read.csv('../../16s_merged/qiime_diversity/otu_table_merged_rarefied/feature-table.tsv', sep='\t', header=T, skip=1, check.names=F)
#arrange columns ABC order
otu_table<-otu_table[ , order(names(otu_table))]
#convert to matrix
otumat<-as.matrix(otu_table[,2:271])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
###-------read in taxonomy table-----------------
tax_table<-read.csv('../../16s_merged/qiime_output/silva_taxonomy_merged.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('../../16s_merged/metadata_16s_merged.csv', header=T,)
#put rows in abc order
sampledata<-sampledata[order(sampledata$id),]
rownames(sampledata)<-sampledata$id
x<-sampledata$id[(which(!(sampledata$id %in% colnames(otumat))))]
sampledata<-sampledata[which(!(sampledata$id %in% x)),] #get rid of samples filtered by qiime (those in x)
sampledata<-sampledata[,2:31] #get rid of id column
###-------read in nwk tree with ape package-----------
tree<-read.tree('../../16s_merged/qiime_output/tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
ps_micro = phyloseq(OTU, TAX, SAM, tree) #note: does not included zero sum otus

#calculate total number of seuqnces
sum(ps_micro@otu_table) #1,085,400

#save ps_micro as file to import and work with microbiome data on its own (this was rarefied by qiime and is good to go!)
saveRDS(ps_micro, file = "micro_physeq.rds")

###----------- next steps------------
#1: import microbiome data -- done!
#2: import diet data -- done!
#filter out non-arthropods from diet dataset
table(ps_diet_raw@tax_table[,2])
chords<-rownames(ps_diet_raw@tax_table[which(!(ps_diet_raw@tax_table[,2]=="p:Arthropoda")),])
non_chords<-rownames(ps_diet_raw@tax_table)
non_chords<-non_chords[which(!(non_chords %in% chords))]
ps_diet_raw2<-prune_taxa(non_chords, ps_diet_raw)
#remove zero sum taxa
ps_diet_raw2<-prune_taxa(taxa_sums(ps_diet_raw2) > 0, ps_diet_raw2)

#3: filter both datasets to include only overlapping samples (for which we have both diet and microbiome data)
####How many samples are in the microbiome dataset?
dim(ps_micro@otu_table) 
####How many samples are in the diet dataset?
dim(ps_diet_raw2@otu_table) 
####There are more samples in the diet dataset. Which ones are also in the microbiome dataset? 
overlap<-colnames(ps_micro@otu_table)[which(colnames(ps_micro@otu_table) %in% colnames(ps_diet_raw@otu_table))] #all micro samples are in diet object
colnames(ps_micro@otu_table)[which(!(colnames(ps_micro@otu_table) %in% colnames(ps_diet_raw@otu_table)))] # none
#4: Remove non-overlapping samples from diet objects
ps_diet<-prune_samples(overlap, ps_diet_raw2)
ps_diet<-prune_taxa(taxa_sums(ps_diet) > 0, ps_diet) #remove zero-sum taxa

#5: rarefy die dataset (microbiome already done)
hist(colSums(ps_micro@otu_table))
summary(colSums(ps_micro@otu_table))
hist(colSums(ps_diet@otu_table))
summary(colSums(ps_diet@otu_table))

rarecurve(t(otu_table(ps_diet)), step=50, lwd=0.8, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,9000))
pdf(file='rarecurve_merged_diet.pdf', height=4, width=4)
rarecurve(t(otu_table(ps_diet)), step=50, lwd=0.7, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,20000), ylim=c(0,200))
dev.off()
pdf(file='figs_diet_and_microbiome_merged/rarecurve_merged_diet2.pdf', height=4, width=4)
rarecurve(t(otu_table(ps_diet)), step=50, lwd=0.7, label=F, xlab='Read count',
          ylab='OTUs',xlim=c(0,20000), ylim=c(0,500))
dev.off()

 ds<-colSums(ps_diet@otu_table)
length(ds[which(ds<15000)])

counts<-ps_micro@sam_data
counts$micro_counts<-ms
counts$diet_counts<-ds

barplot(table(counts$Species), las=2, main='original', ylim=c(0,25))
barplot(table(counts$Species[which(counts$diet_counts>15000)]), las=2, main='15000', ylim=c(0,25))
abline(h=10, col='red', lty=2) #most still have >10 indiv/species


#5.1 create a new object that holds all of the samples below both thresholds. We will use this to filter them from the diet and microbiome phyloseq objects.
####using the counts dataframe, find the 54 samples below diet threshold of 15000. Assign them to a new object called t1. 
t1<-rownames(counts)[which(counts$diet_counts<15000)]
###We will use t3 to prune our two phyloseq objects, and then rarefy the diet one!

#prune low read count diet samples
micro_rare<-prune_samples(rownames(counts)[which(!(rownames(counts) %in% t1))], ps_micro) #already rarefied, name to match diet below
ps_diet2<-prune_samples(rownames(counts)[which(!(rownames(counts) %in% t1))], ps_diet)

#rarefy diet dataset (default is lowest read count, which is 4020 and 15118)
summary(colSums(ps_diet2@otu_table))
diet_rare<-rarefy_even_depth(ps_diet2, rngseed=8, replace=F)

#check that non-zero sum taxa removed
table(taxa_sums(diet_rare) > 0) #check!
table(taxa_sums(micro_rare) > 0) #need to remove
micro_rare<-prune_taxa(taxa_sums(micro_rare) > 0, micro_rare) #remove zero-sum taxa

#root the diet tree!
arachnids<-rownames(diet_rare@tax_table)[which(diet_rare@tax_table[,3]=='c:Arachnida')] #build vector of arachnid OTUs
sample(arachnids, 1) #choose random arachnid to root the tree on, got OTU1598
diet_rare@phy_tree<-root(diet_rare@phy_tree, outgroup = 'OTU1598', resolve.root=T)

#6: calculate alpha and beta for each dataset
alpha_micro<-diversity(micro_rare@otu_table, index='shannon', MARGIN=2)
alpha_diet<-diversity(diet_rare@otu_table, index='shannon', MARGIN=2)
plot(alpha_diet, alpha_micro)
abline(lm(alpha_micro~alpha_diet))

div<-micro_rare@sam_data
div$micro_shannon<-alpha_micro
div$diet_shannon<-alpha_diet

library(randomcoloR)
palette(distinctColorPalette(15))
par(mar=c(5, 4, 2, 7.3), xpd=TRUE)
plot(div$micro_shannon, div$diet_shannon, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (Shannon)', ylab='Diet alpha diversity (Shannon)')
abline(lm(div$diet_shannon~div$micro_shannon), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))


plot(div$micro_shannon[which(div$Species=='OVEN')], div$diet_shannon[which(div$Species=='OVEN')])
plot(div$micro_shannon[which(div$Species=='BTNW')], div$diet_shannon[which(div$Species=='BTNW')])
plot(div$micro_shannon[which(div$Species=='BTBW')], div$diet_shannon[which(div$Species=='BTBW')])
plot(div$micro_shannon[which(div$Species=='BAWW')], div$diet_shannon[which(div$Species=='BAWW')])

alpha_micro<-estimate_richness(micro_rare, split=T, measures=c('Observed', 'Shannon'))
alpha_diet<-estimate_richness(diet_rare, split=T, measures=c('Observed', 'Shannon'))

plot(alpha_diet$Observed, alpha_micro$Observed)
abline(lm(alpha_micro$Observed~alpha_diet$Observed))

div$observed_micro<-alpha_micro$Observed
div$observed_diet<-alpha_diet$Observed

#7: analyze diet and microbiome diversity together to test our predictions!
plot_richness(diet_rare, x="Species", color="Species",measures="Observed") + theme(legend.position="none") +
  geom_boxplot() 

#individual species plots alpha diet~alpha microbiome
plot(div$micro_shannon[which(div$Species=='AMRE')], div$diet_shannon[which(div$Species=='AMRE')], main='AMRE', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='AMRE')]~div$micro_shannon[which(div$Species=='AMRE')]))

plot(div$micro_shannon[which(div$Species=='BAWW')], div$diet_shannon[which(div$Species=='BAWW')], main='BAWW', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='BAWW')]~div$micro_shannon[which(div$Species=='BAWW')]))

plot(div$micro_shannon[which(div$Species=='BLBW')], div$diet_shannon[which(div$Species=='BLBW')], main='BLBW', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='BLBW')]~div$micro_shannon[which(div$Species=='BLBW')]))

plot(div$micro_shannon[which(div$Species=='BTBW')], div$diet_shannon[which(div$Species=='BTBW')], main='BTBW', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='BTBW')]~div$micro_shannon[which(div$Species=='BTBW')]))

plot(div$micro_shannon[which(div$Species=='BTNW')], div$diet_shannon[which(div$Species=='BTNW')], main='BTNW', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='BTNW')]~div$micro_shannon[which(div$Species=='BTNW')]))

plot(div$micro_shannon[which(div$Species=='CAWA')], div$diet_shannon[which(div$Species=='CAWA')], main='CAWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='CAWA')]~div$micro_shannon[which(div$Species=='CAWA')]))

plot(div$micro_shannon[which(div$Species=='COYE')], div$diet_shannon[which(div$Species=='COYE')], main='COYE', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='COYE')]~div$micro_shannon[which(div$Species=='COYE')]))

plot(div$micro_shannon[which(div$Species=='CSWA')], div$diet_shannon[which(div$Species=='CSWA')], main='CSWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='CSWA')]~div$micro_shannon[which(div$Species=='CSWA')]))

plot(div$micro_shannon[which(div$Species=='HOWA')], div$diet_shannon[which(div$Species=='HOWA')], main='HOWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='HOWA')]~div$micro_shannon[which(div$Species=='HOWA')]))

plot(div$micro_shannon[which(div$Species=='MAWA')], div$diet_shannon[which(div$Species=='MAWA')], main='MAWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='MAWA')]~div$micro_shannon[which(div$Species=='MAWA')]))

plot(div$micro_shannon[which(div$Species=='MYWA')], div$diet_shannon[which(div$Species=='MYWA')], main='MYWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='MYWA')]~div$micro_shannon[which(div$Species=='MYWA')]))

plot(div$micro_shannon[which(div$Species=='NAWA')], div$diet_shannon[which(div$Species=='NAWA')], main='NAWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='NAWA')]~div$micro_shannon[which(div$Species=='NAWA')]))

plot(div$micro_shannon[which(div$Species=='NOPA')], div$diet_shannon[which(div$Species=='NOPA')], main='NOPA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='NOPA')]~div$micro_shannon[which(div$Species=='NOPA')]))

plot(div$micro_shannon[which(div$Species=='OVEN')], div$diet_shannon[which(div$Species=='OVEN')], main='OVEN', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='OVEN')]~div$micro_shannon[which(div$Species=='OVEN')]))

plot(div$micro_shannon[which(div$Species=='WEWA')], div$diet_shannon[which(div$Species=='WEWA')], main='WEWA', xlab='Microbiome Alpha Diversity', ylab='Diet Alpha Diversity')
abline(lm(div$diet_shannon[which(div$Species=='WEWA')]~div$micro_shannon[which(div$Species=='WEWA')]))

###--find mean alpha diversity for each species, diet and micro
a1<-mean(div$micro_shannon[which(div$Species=='AMRE')])
b1<-mean(div$diet_shannon[which(div$Species=='AMRE')])
c1<-mean(div$micro_chao1[which(div$Species=='AMRE')])
d1<-mean(div$diet_chao1[which(div$Species=='AMRE')])
e1<-mean(div$micro_PD[which(div$Species=='AMRE')])
f1<-mean(div$diet_PD[which(div$Species=='AMRE')])
#95% CIs
ca1<-CI(div$micro_shannon[which(div$Species=='AMRE')], ci=0.95)
cd1<-CI(div$diet_shannon[which(div$Species=='AMRE')], ci=0.95)

a2<-mean(div$micro_shannon[which(div$Species=='BAWW')])
b2<-mean(div$diet_shannon[which(div$Species=='BAWW')])
c2<-mean(div$micro_chao1[which(div$Species=='BAWW')])
d2<-mean(div$diet_chao1[which(div$Species=='BAWW')])
e2<-mean(div$micro_PD[which(div$Species=='BAWW')])
f2<-mean(div$diet_PD[which(div$Species=='BAWW')])
ca2<-CI(div$micro_shannon[which(div$Species=='BAWW')], ci=0.95)
cd2<-CI(div$diet_shannon[which(div$Species=='BAWW')], ci=0.95)

a3<-mean(div$micro_shannon[which(div$Species=='BLBW')])
b3<-mean(div$diet_shannon[which(div$Species=='BLBW')])
c3<-mean(div$micro_chao1[which(div$Species=='BLBW')])
d3<-mean(div$diet_chao1[which(div$Species=='BLBW')])
e3<-mean(div$micro_PD[which(div$Species=='BLBW')])
f3<-mean(div$diet_PD[which(div$Species=='BLBW')])
ca3<-CI(div$micro_shannon[which(div$Species=='BLBW')], ci=0.95)
cd3<-CI(div$diet_shannon[which(div$Species=='BLBW')], ci=0.95)

a4<-mean(div$micro_shannon[which(div$Species=='BTBW')])
b4<-mean(div$diet_shannon[which(div$Species=='BTBW')])
c4<-mean(div$micro_chao1[which(div$Species=='BTBW')])
d4<-mean(div$diet_chao1[which(div$Species=='BTBW')])
e4<-mean(div$micro_PD[which(div$Species=='BTBW')])
f4<-mean(div$diet_PD[which(div$Species=='BTBW')])
ca4<-CI(div$micro_shannon[which(div$Species=='BTBW')], ci=0.95)
cd4<-CI(div$diet_shannon[which(div$Species=='BTBW')], ci=0.95)

a5<-mean(div$micro_shannon[which(div$Species=='BTNW')])
b5<-mean(div$diet_shannon[which(div$Species=='BTNW')])
c5<-mean(div$micro_chao1[which(div$Species=='BTNW')])
d5<-mean(div$diet_chao1[which(div$Species=='BTNW')])
e5<-mean(div$micro_PD[which(div$Species=='BTNW')])
f5<-mean(div$diet_PD[which(div$Species=='BTNW')])
ca5<-CI(div$micro_shannon[which(div$Species=='BTNW')], ci=0.95)
cd5<-CI(div$diet_shannon[which(div$Species=='BTNW')], ci=0.95)

a6<-mean(div$micro_shannon[which(div$Species=='CAWA')])
b6<-mean(div$diet_shannon[which(div$Species=='CAWA')])
c6<-mean(div$micro_chao1[which(div$Species=='CAWA')])
d6<-mean(div$diet_chao1[which(div$Species=='CAWA')])
e6<-mean(div$micro_PD[which(div$Species=='CAWA')])
f6<-mean(div$diet_PD[which(div$Species=='CAWA')])
ca6<-CI(div$micro_shannon[which(div$Species=='CAWA')], ci=0.95)
cd6<-CI(div$diet_shannon[which(div$Species=='CAWA')], ci=0.95)

a7<-mean(div$micro_shannon[which(div$Species=='COYE')])
b7<-mean(div$diet_shannon[which(div$Species=='COYE')])
c7<-mean(div$micro_chao1[which(div$Species=='COYE')])
d7<-mean(div$diet_chao1[which(div$Species=='COYE')])
e7<-mean(div$micro_PD[which(div$Species=='COYE')])
f7<-mean(div$diet_PD[which(div$Species=='COYE')])
ca7<-CI(div$micro_shannon[which(div$Species=='COYE')], ci=0.95)
cd7<-CI(div$diet_shannon[which(div$Species=='COYE')], ci=0.95)

a8<-mean(div$micro_shannon[which(div$Species=='CSWA')])
b8<-mean(div$diet_shannon[which(div$Species=='CSWA')])
c8<-mean(div$micro_chao1[which(div$Species=='CSWA')])
d8<-mean(div$diet_chao1[which(div$Species=='CSWA')])
e8<-mean(div$micro_PD[which(div$Species=='CSWA')])
f8<-mean(div$diet_PD[which(div$Species=='CSWA')])
ca8<-CI(div$micro_shannon[which(div$Species=='CSWA')], ci=0.95)
cd8<-CI(div$diet_shannon[which(div$Species=='CSWA')], ci=0.95)

a9<-mean(div$micro_shannon[which(div$Species=='HOWA')])
b9<-mean(div$diet_shannon[which(div$Species=='HOWA')])
c9<-mean(div$micro_chao1[which(div$Species=='HOWA')])
d9<-mean(div$diet_chao1[which(div$Species=='HOWA')])
e9<-mean(div$micro_PD[which(div$Species=='HOWA')])
f9<-mean(div$diet_PD[which(div$Species=='HOWA')])
ca9<-CI(div$micro_shannon[which(div$Species=='HOWA')], ci=0.95)
cd9<-CI(div$diet_shannon[which(div$Species=='HOWA')], ci=0.95)

a10<-mean(div$micro_shannon[which(div$Species=='MAWA')])
b10<-mean(div$diet_shannon[which(div$Species=='MAWA')])
c10<-mean(div$micro_chao1[which(div$Species=='MAWA')])
d10<-mean(div$diet_chao1[which(div$Species=='MAWA')])
e10<-mean(div$micro_PD[which(div$Species=='MAWA')])
f10<-mean(div$diet_PD[which(div$Species=='MAWA')])
ca10<-CI(div$micro_shannon[which(div$Species=='MAWA')], ci=0.95)
cd10<-CI(div$diet_shannon[which(div$Species=='MAWA')], ci=0.95)

a11<-mean(div$micro_shannon[which(div$Species=='MYWA')])
b11<-mean(div$diet_shannon[which(div$Species=='MYWA')])
c11<-mean(div$micro_chao1[which(div$Species=='MYWA')])
d11<-mean(div$diet_chao1[which(div$Species=='MYWA')])
e11<-mean(div$micro_PD[which(div$Species=='MYWA')])
f11<-mean(div$diet_PD[which(div$Species=='MYWA')])
ca11<-CI(div$micro_shannon[which(div$Species=='MYWA')], ci=0.95)
cd11<-CI(div$diet_shannon[which(div$Species=='MYWA')], ci=0.95)

a12<-mean(div$micro_shannon[which(div$Species=='NAWA')])
b12<-mean(div$diet_shannon[which(div$Species=='NAWA')])
c12<-mean(div$micro_chao1[which(div$Species=='NAWA')])
d12<-mean(div$diet_chao1[which(div$Species=='NAWA')])
e12<-mean(div$micro_PD[which(div$Species=='NAWA')])
f12<-mean(div$diet_PD[which(div$Species=='NAWA')])
ca12<-CI(div$micro_shannon[which(div$Species=='NAWA')], ci=0.95)
cd12<-CI(div$diet_shannon[which(div$Species=='NAWA')], ci=0.95)

a13<-mean(div$micro_shannon[which(div$Species=='NOPA')])
b13<-mean(div$diet_shannon[which(div$Species=='NOPA')])
c13<-mean(div$micro_chao1[which(div$Species=='NOPA')])
d13<-mean(div$diet_chao1[which(div$Species=='NOPA')])
e13<-mean(div$micro_PD[which(div$Species=='NOPA')])
f13<-mean(div$diet_PD[which(div$Species=='NOPA')])
ca13<-CI(div$micro_shannon[which(div$Species=='NOPA')], ci=0.95)
cd13<-CI(div$diet_shannon[which(div$Species=='NOPA')], ci=0.95)

a14<-mean(div$micro_shannon[which(div$Species=='OVEN')])
b14<-mean(div$diet_shannon[which(div$Species=='OVEN')])
c14<-mean(div$micro_chao1[which(div$Species=='OVEN')])
d14<-mean(div$diet_chao1[which(div$Species=='OVEN')])
e14<-mean(div$micro_PD[which(div$Species=='OVEN')])
f14<-mean(div$diet_PD[which(div$Species=='OVEN')])
ca14<-CI(div$micro_shannon[which(div$Species=='OVEN')], ci=0.95)
cd14<-CI(div$diet_shannon[which(div$Species=='OVEN')], ci=0.95)

a15<-mean(div$micro_shannon[which(div$Species=='WEWA')])
b15<-mean(div$diet_shannon[which(div$Species=='WEWA')])
c15<-mean(div$micro_chao1[which(div$Species=='WEWA')])
d15<-mean(div$diet_chao1[which(div$Species=='WEWA')])
e15<-mean(div$micro_PD[which(div$Species=='WEWA')])
f15<-mean(div$diet_PD[which(div$Species=='WEWA')])
ca15<-CI(div$micro_shannon[which(div$Species=='WEWA')], ci=0.95)
cd15<-CI(div$diet_shannon[which(div$Species=='WEWA')], ci=0.95)

#calculate beta diversity!
bray1a<-mean(distance(subset_samples(micro_rare, Species=='AMRE'), "bray"))
bray1b<-mean(distance(subset_samples(diet_rare, Species=='AMRE'), "bray"))
bray2a<-mean(distance(subset_samples(micro_rare, Species=='BAWW'), "bray"))
bray2b<-mean(distance(subset_samples(diet_rare, Species=='BAWW'), "bray"))
bray3a<-mean(distance(subset_samples(micro_rare, Species=='BLBW'), "bray"))
bray3b<-mean(distance(subset_samples(diet_rare, Species=='BLBW'), "bray"))
bray4a<-mean(distance(subset_samples(micro_rare, Species=='BTBW'), "bray"))
bray4b<-mean(distance(subset_samples(diet_rare, Species=='BTBW'), "bray"))
bray5a<-mean(distance(subset_samples(micro_rare, Species=='BTNW'), "bray"))
bray5b<-mean(distance(subset_samples(diet_rare, Species=='BTNW'), "bray"))
bray6a<-mean(distance(subset_samples(micro_rare, Species=='CAWA'), "bray"))
bray6b<-mean(distance(subset_samples(diet_rare, Species=='CAWA'), "bray"))
bray7a<-mean(distance(subset_samples(micro_rare, Species=='COYE'), "bray"))
bray7b<-mean(distance(subset_samples(diet_rare, Species=='COYE'), "bray"))
bray8a<-mean(distance(subset_samples(micro_rare, Species=='CSWA'), "bray"))
bray8b<-mean(distance(subset_samples(diet_rare, Species=='CSWA'), "bray"))
bray9a<-mean(distance(subset_samples(micro_rare, Species=='HOWA'), "bray"))
bray9b<-mean(distance(subset_samples(diet_rare, Species=='HOWA'), "bray"))
bray10a<-mean(distance(subset_samples(micro_rare, Species=='MAWA'), "bray"))
bray10b<-mean(distance(subset_samples(diet_rare, Species=='MAWA'), "bray"))
bray11a<-mean(distance(subset_samples(micro_rare, Species=='MYWA'), "bray"))
bray11b<-mean(distance(subset_samples(diet_rare, Species=='MYWA'), "bray"))
bray12a<-mean(distance(subset_samples(micro_rare, Species=='NAWA'), "bray"))
bray12b<-mean(distance(subset_samples(diet_rare, Species=='NAWA'), "bray"))
bray13a<-mean(distance(subset_samples(micro_rare, Species=='NOPA'), "bray"))
bray13b<-mean(distance(subset_samples(diet_rare, Species=='NOPA'), "bray"))
bray14a<-mean(distance(subset_samples(micro_rare, Species=='OVEN'), "bray"))
bray14b<-mean(distance(subset_samples(diet_rare, Species=='OVEN'), "bray"))
bray15a<-mean(distance(subset_samples(micro_rare, Species=='WEWA'), "bray"))
bray15b<-mean(distance(subset_samples(diet_rare, Species=='WEWA'), "bray"))


bray_micro<-distance(micro_rare, 'bray')
bray_diet<-distance(diet_rare, 'bray')

#unifrac
unifrac1a<-mean(distance(subset_samples(micro_rare,Species=='AMRE'),'unifrac'))
unifrac1b<-mean(distance(subset_samples(diet_rare,Species=='AMRE'),'unifrac'))

unifrac2a<-mean(distance(subset_samples(micro_rare,Species=='BAWW'),'unifrac'))
unifrac2b<-mean(distance(subset_samples(diet_rare,Species=='BAWW'),'unifrac'))

unifrac3a<-mean(distance(subset_samples(micro_rare,Species=='BLBW'),'unifrac'))
unifrac3b<-mean(distance(subset_samples(diet_rare,Species=='BLBW'),'unifrac'))

unifrac4a<-mean(distance(subset_samples(micro_rare,Species=='BTBW'),'unifrac'))
unifrac4b<-mean(distance(subset_samples(diet_rare,Species=='BTBW'),'unifrac'))

unifrac5a<-mean(distance(subset_samples(micro_rare,Species=='BTNW'),'unifrac'))
unifrac5b<-mean(distance(subset_samples(diet_rare,Species=='BTNW'),'unifrac'))

unifrac6a<-mean(distance(subset_samples(micro_rare,Species=='CAWA'),'unifrac'))
unifrac6b<-mean(distance(subset_samples(diet_rare,Species=='CAWA'),'unifrac'))

unifrac7a<-mean(distance(subset_samples(micro_rare,Species=='COYE'),'unifrac'))
unifrac7b<-mean(distance(subset_samples(diet_rare,Species=='COYE'),'unifrac'))

unifrac8a<-mean(distance(subset_samples(micro_rare,Species=='CSWA'),'unifrac'))
unifrac8b<-mean(distance(subset_samples(diet_rare,Species=='CSWA'),'unifrac'))

unifrac9a<-mean(distance(subset_samples(micro_rare,Species=='HOWA'),'unifrac'))
unifrac9b<-mean(distance(subset_samples(diet_rare,Species=='HOWA'),'unifrac'))

unifrac10a<-mean(distance(subset_samples(micro_rare,Species=='MAWA'),'unifrac'))
unifrac10b<-mean(distance(subset_samples(diet_rare,Species=='MAWA'),'unifrac'))

unifrac11a<-mean(distance(subset_samples(micro_rare,Species=='MYWA'),'unifrac'))
unifrac11b<-mean(distance(subset_samples(diet_rare,Species=='MYWA'),'unifrac'))

unifrac12a<-mean(distance(subset_samples(micro_rare,Species=='NAWA'),'unifrac'))
unifrac12b<-mean(distance(subset_samples(diet_rare,Species=='NAWA'),'unifrac'))

unifrac13a<-mean(distance(subset_samples(micro_rare,Species=='NOPA'), 'unifrac'))
unifrac13b<-mean(distance(subset_samples(diet_rare,Species=='NOPA'), 'unifrac'))

unifrac14a<-mean(distance(subset_samples(micro_rare,Species=='OVEN'),'unifrac'))
unifrac14b<-mean(distance(subset_samples(diet_rare,Species=='OVEN'),'unifrac'))

unifrac15a<-mean(distance(subset_samples(micro_rare,Species=='WEWA'),'unifrac'))
unifrac15b<-mean(distance(subset_samples(diet_rare,Species=='WEWA'),'unifrac'))


#build dataframe for diet index score!
d_up<-CI(c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15), ci=0.95)

index_ab<-data.frame("alpha_diet"=c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15),
                     "alpha_diet_upper"=c(cd1[1],cd2[1],cd3[1],cd4[1],cd5[1],cd6[1],cd7[1],cd8[1],cd9[1],cd10[1],cd11[1],cd12[1],cd13[1],cd14[1],cd15[1]),
                     "alpha_diet_lower"=c(cd1[3],cd2[3],cd3[3],cd4[3],cd5[3],cd6[3],cd7[3],cd8[3],cd9[3],cd10[3],cd11[3],cd12[3],cd13[3],cd14[3],cd15[3]),
                     "alpha_micro"=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15),
                     "alpha_micro_upper"=c(ca1[1],ca2[1],ca3[1],ca4[1],ca5[1],ca6[1],ca7[1],ca8[1],ca9[1],ca10[1],ca11[1],ca12[1],ca13[1],ca14[1],ca15[1]),
                     "alpha_micro_lower"=c(ca1[3],ca2[3],ca3[3],ca4[3],ca5[3],ca6[3],ca7[3],ca8[3],ca9[3],ca10[3],ca11[3],ca12[3],ca13[3],ca14[3],ca15[3]),
                     "beta_diet"=c(bray1b,bray2b,bray3b,bray4b,bray5b,bray6b,bray7b,bray8b,bray9b,bray10b,bray11b,bray12b,bray13b,bray14b,bray15b),
                     "beta_micro"=c(bray1a,bray2a,bray3a,bray4a,bray5a,bray6a,bray7a,bray8a,bray9a,bray10a,bray11a,bray12a,bray13a,bray14a,bray15a),
                     'unifrac_diet'=c(unifrac1b,unifrac2b,unifrac3b,unifrac4b,unifrac5b,unifrac6b,unifrac7b,unifrac8b,unifrac9b,unifrac10b,unifrac11b,unifrac12b,unifrac13b,unifrac14b,unifrac15b),
                     'unifrac_micro'=c(unifrac1a,unifrac2a,unifrac3a,unifrac4a,unifrac5a,unifrac6a,unifrac7a,unifrac8a,unifrac9a,unifrac10a,unifrac11a,unifrac12a,unifrac13a,unifrac14a,unifrac15a))
rownames(index_ab)<-c('AMRE','BAWW','BLBW','BTBW','BTNW','CAWA','COYE','CSWA','HOWA','MAWA','MYWA','NAWA','NOPA','OVEN','WEWA')

plot(index_ab$alpha_diet, index_ab$beta_diet)
plot(index_ab$alpha_diet, index_ab$unifrac_diet)
plot(index_ab$alpha_diet, index_ab$alpha_micro)
abline(lm(index_ab$alpha_micro~index_ab$alpha_diet))

summary(lm(index_ab$alpha_micro~index_ab$alpha_diet))

index_ab$i1<-index_ab$alpha_diet+index_ab$beta_diet
index_ab$i2<-index_ab$alpha_diet+index_ab$unifrac_diet

plot(index_ab$i2, index_ab$alpha_micro)
abline(lm(index_ab$alpha_micro~index_ab$i2))

hist(index_ab$i1)
hist(index_ab$i2)
plot(index_ab$i1~index_ab$i2)

pdf(file='figs_diet_and_microbiome_merged/hist_index_scores2017-2020.pdf', height=6, width=5)
par(mfrow=c(2,1), mar=c(4,4,2,3))
hist(index_ab$i1, breaks=8, xlab='Index score', main='Index of diet specialization using BC distance')
hist(index_ab$i2, breaks=8, xlab='Index score', main='Index of diet specialization using UniFrac distance')
dev.off()

pdf(file='figs_diet_and_microbiome_merged/plot_index_scores2017-2020.pdf', height=5, width=4)
plot(index_ab$i1~index_ab$i2, ylab='Diet index using Bray-Curtis distance',
     xlab='Diet index using UniFrac distance', main='2017-2020')
dev.off()

index_ab[order(index_ab$i1),]

diet_rare@sam_data$Sequencing_lane<-as.factor(diet_rare@sam_data$Sequencing_lane) #make non-continuous
ordinate(diet_rare, "PCoA", "bray") %>% 
  plot_ordination(diet_rare, ., color='Sequencing_lane', title = "Diet Bray-Curtis", type='samples') +
  geom_point(size=2)
ordinate(diet_rare, "PCoA", "bray") %>% 
  plot_ordination(diet_rare, ., color='Species', title = "Diet Bray-Curtis", type='samples') +
  geom_point(size=2)


ps0<-tax_glom(diet_rare, taxrank = "Class", NArm=F) #glom by class
mps0<-psmelt(ps0)
table(mps0$Class)
ggplot(mps0, aes(x = sample_Species, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('cornsilk3','darkblue',
                             'cornflowerblue','coral3','darkgreen',
                             'darkolivegreen3','black')) +
  ggtitle('Diet') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))

###diet relative abundance plot distinguishing only top 7 orders
topc<-as.data.frame(diet_rare@tax_table)
sort(table(topc$Order), decreasing=T)
topcc<-names(sort(table(topc$Order), decreasing=T)[1:7]) #top 7 orders

ps1<-tax_glom(diet_rare, taxrank = "Order", NArm=F) #glom by order
mps1<-psmelt(ps1)
mps1$Order[which(!(mps1$Order %in% topcc))]<-'other' #change non-common orders to other 
table(mps1$Order)

mps1$Order2<-as.factor(mps1$Order)
mps1$Order2<-relevel(mps1$Order2, 'other')

mps1$Order[which(mps1$Class=='c:Arachnida' & mps1$Abundance>0)]

pdf(file='figs_diet_and_microbiome_merged/abundance_diet_2017-2020.pdf', height=5, width=4)
level_order<-c('BTNW','MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','HOWA','CAWA','COYE','NAWA','BAWW','WEWA','OVEN')
ggplot(mps1, aes(x = factor(sample_Species, level=level_order), y = Abundance, fill = Order2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','aquamarine3',
                             'aquamarine4','darkseagreen','palegreen4',
                             'wheat','tan','dimgray' )) +
  ggtitle('') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(axis.text = element_text(size = 11)) +
  theme(legend.position="none")
dev.off()

###microbiome relative abundance top 8 phyla--don't use, use dataset that didn't filter out 15K diet
topm<-as.data.frame(micro_rare@tax_table)
sort(table(topm$Phylum), decreasing=T)
topmm<-names(sort(table(topm$Phylum), decreasing=T)[1:8]) #top 8 orders

ps2<-tax_glom(micro_rare, taxrank = "Phylum", NArm=F) #glom by phylum
mps2<-psmelt(ps2)
mps2$Phylum[which(!(mps2$Phylum %in% topmm))]<-'other' #change non-common phyla to other 
table(mps2$Phylum)

mps2$Phylum2<-as.factor(mps2$Phylum)
mps2$Phylum2<-relevel(mps2$Phylum2, 'other')
ggplot(mps2, aes(x = sample_Species, y = Abundance, fill = Phylum2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('cornsilk3','darkgray',
                             'coral1','orange','brown','darkred',
                             'red','pink','black')) +
  ggtitle('Microbiome') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0))



#calculate total number of seuqnces
sum(ps_diet_raw@otu_table) #17,869,484
sum(diet_rare@otu_table) #3,265,488 
sum(ps_micro@otu_table) #1,085,400 
sum(micro_rare@otu_table) #868,320 

#KW test alpha by species
kruskal.test(alpha_diet$Observed, sample_data(diet_rare)$Species)
kruskal.test(alpha_micro$Observed, sample_data(micro_rare)$Species)
kruskal.test(alpha_diet$Shannon, sample_data(diet_rare)$Species)
kruskal.test(alpha_micro$Shannon, sample_data(micro_rare)$Species)

#PERMANOVAs
#generate new distance matrices
#---don't need to do since singleton spp removed
#bray_micro2<-distance(micro_rare2, 'bray')
#bray_diet2<-distance(diet_rare2, 'bray')

md<-data.frame(micro_rare@sam_data, row.names=rownames(micro_rare@sam_data))
md$Year<-as.factor(md$Year) #make non-continuous
md$Sequencing_lane<-as.factor(md$Sequencing_lane) #make non-continuous

adonis(bray_micro ~ Species, data=md) #ns
adonis(bray_micro ~ Year, data=md)
adonis(bray_micro ~ State, data=md)
adonis(bray_micro ~ Sequencing_lane, data=md) #r2=0.15775, p=0.001***
adonis(bray_micro ~ Diet_type, data=md) #r2=0.00917,  p=0.435
adonis(bray_micro ~ Habitat, data=md) #r2=0.0197  p=0.279

unifrac_micro<-distance(micro_rare,'unifrac')
unifrac_diet<-distance(diet_rare,'unifrac')

adonis(unifrac_micro~Species,data=md)
adonis(unifrac_micro~Year,data=md)
adonis(unifrac_micro~State,data=md)
adonis(unifrac_micro~Sequencing_lane,data=md) #r2=0.09061, p=0.001
adonis(unifrac_micro~Diet_type, data=md) #r2=0.01003,  p=0.191
adonis(unifrac_micro~Habitat, data=md) #ns

wunifrac_micro<-distance(micro_rare,'wunifrac')
adonis(wunifrac_micro~Diet_type, data=md) #r2=0.01341,  p=0.19

index_ab
hist(index_ab$i1)
index_ab$Diet_Type<-c('Specialist','Generalist','Intermediate',
                      'Intermediate','Intermediate','Generalist',
                      'Generalist','Specialist','Generalist',
                      'Intermediate','Intermediate','Intermediate',
                      'Intermediate','Intermediate','Specialist')

md$Diet_type<-index_ab$Diet_Type[match(md$Species, rownames(index_ab))]

#broad assumption of habitat type...perhaps too broad
index_ab$Habitat<-c('Mixed','Forest, Mixed','Forest, High Canopy',
                    'Forest, Low Canopy','Forest, High Canopy',
                    'Forest, Low Canopy','Open','Open','Forest, Mixed',
                    'Forest, Mixed','Forest, Mixed','Mixed','Open','Forest, Low Canopy', 'Forest, Low Canopy')
#save index_ab table to be read in am map habitat/diet type
saveRDS(index_ab, file = "index_ab.rds")

md$Habitat<-index_ab$Habitat[match(md$Species, rownames(index_ab))]


micro_rare@sam_data$Diet_type<-md$Diet_type
micro_rare@sam_data$Habitat<-md$Habitat

ordinate(micro_rare, "PCoA", "bray") %>% 
  plot_ordination(micro_rare, ., color='Diet_type', title = "Microbiome Bray-Curtis", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Diet type"))
ordinate(micro_rare, "PCoA", "unifrac") %>% 
  plot_ordination(micro_rare, ., color='Diet_type', title = "Micrbiome UniFrac", type='samples') +
  geom_point(size=2)
ordinate(micro_rare, "PCoA", "bray") %>% 
  plot_ordination(micro_rare, ., color='Habitat', title = "Microbiome Bray-Curtis", type='samples') +
  geom_point(size=2)
ordinate(micro_rare, "PCoA", "unifrac") %>% 
  plot_ordination(micro_rare, ., color='Habitat', title = "Microbiome Bray-Curtis", type='samples') +
  geom_point(size=2)

#pdf(file='figs_diet_and_microbiome_merged/boxplot_shannon_diettype2017-2020.pdf', height=4, width=3.5)
plot_richness(micro_rare, x="Diet_type", color="Diet_type",measures="Shannon") + theme(legend.position="none") +
  geom_boxplot(alpha=0.7, width=0.5) +scale_color_manual(values=c('#CC79A7', '#E69F00', '#56B4E9')) +
  theme_bw() +
  xlab('Diet type') + ylab('Microbiome shannon diversity') +
  theme(axis.text = element_text(size = 11)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(legend.position="none")
dev.off()

plot_richness(micro_rare, x="Diet_type", color="Diet_type",measures="Shannon") + theme(legend.position="none") +
  geom_boxplot(alpha=0.7, width=0.5) +scale_color_manual(values=c('chocolate4', 'aquamarine4', 'darkgoldenrod')) +
  theme_bw()

kruskal.test(alpha_micro$Shannon, sample_data(micro_rare)$Diet_type) #ns chi=0.013663, p=0.9932
kruskal.test(alpha_micro$Observed, sample_data(micro_rare)$Diet_type) #ns

extremes<-as.data.frame(div)
extremes$Diet_type<-md$Diet_type
extremes<-extremes[which(!(extremes$Diet_type=='Intermediate')),] #93 individuals, just generalists and specailists
kruskal.test(extremes$micro_shannon, extremes$Diet_type) #ns chi-squared=0.048774, df=1, p=0.8252, not significant even if intermediates are excluded


md$shannon<-alpha_micro
boxplot(md$shannon$Shannon[which(md$Diet_type=='Generalist')], md$shannon$Shannon[which(md$Diet_type=='Intermediate')], md$shannon$Shannon[which(md$Diet_type=='Specialist')],
        boxwex=0.5, names=c("Generalist",'Intermediate',"Specialist"),ylab='Shannon index (alpha)', outline=F)
stripchart(md$shannon$Shannon[which(md$Diet_type=='Generalist')],vertical=T,pch=21, method='jitter',add=T,bg='blue')
stripchart(md$shannon$Shannon[which(md$Diet_type=='Intermediate')],vertical=T,pch=21, method='jitter',add=T,bg='red', at=2)
stripchart(md$shannon$Shannon[which(md$Diet_type=='Specialist')],vertical=T,pch=21, method='jitter',add=T,bg='yellow', at=3)


library(randomcoloR)
palette(distinctColorPalette(15))
#pdf(file='~/Documents/Conferences:Presentations/AOS-SCO2021/diet_micro_scatter.pdf', height=5, width=5)
pdf(file='figs_diet_and_microbiome_merged/diet_micro_scatter2017-2020.pdf', height=5, width=5)
par(mar=c(5, 4, 2, 6.5), xpd=TRUE)
plot(div$micro_shannon, div$diet_shannon, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (Shannon)', ylab='Diet alpha diversity (Shannon)',
     main='2017-2020')
abline(lm(div$diet_shannon~div$micro_shannon), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(div$micro_shannon, div$diet_shannon, method='kendall') #tau=0.005426357, p=0.9055

#mean alpha by species
#pdf(file='~/Documents/Conferences:Presentations/AOS-SCO2021/diet_micro_scatter_means.pdf', height=5, width=5)
par(mar=c(5, 4, 2, 6.5), xpd=TRUE)
plot(index_ab$alpha_micro, index_ab$alpha_diet, col=factor(rownames(index_ab)), lwd=2,
     xlab='Mean microbiome alpha diversity (Shannon)', ylab='Mean diet alpha diversity (Shannon)')
abline(lm(index_ab$alpha_diet~index_ab$alpha_micro), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(rownames(index_ab))), col=1:nlevels(factor(rownames(index_ab))), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(index_ab$alpha_micro, index_ab$alpha_diet, method='kendall')

ggplot(index_ab, aes(alpha_micro,alpha_diet)) +
  geom_point() +
  geom_errorbar(aes(ymin=alpha_diet_lower, ymax=alpha_diet_upper),linetype=2, size=0.3) +
  geom_errorbarh(aes(xmin=alpha_micro_lower, xmax=alpha_micro_upper), linetype=2, size=0.3)


#alpha micro by diet index
plot(index_ab$alpha_micro, index_ab$i2, col=factor(rownames(index_ab)), lwd=2,
     xlab='Mean microbiome alpha diversity (Shannon)', ylab='Index of diet specialization')
abline(lm(index_ab$i2~index_ab$alpha_micro), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(rownames(index_ab))), col=1:nlevels(factor(rownames(index_ab))), pch=1, pt.lwd=2, inset=c(-0.35,0))
cor.test(index_ab$alpha_micro, index_ab$i2, method='kendall')



library(metagMisc)
library(PhyloMeasures)
div$micro_pd<-phyloseq_phylo_div(micro_rare, measures='PD')
kruskal.test(div$micro_pd$PD, div$Species)

#pdf(file='~/Documents/Conferences:Presentations/AOS-SCO2021/abundance_diet.pdf', height=5, width=5)
ggplot(mps1, aes(x = sample_Species, y = Abundance, fill = Order2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','aquamarine3',
                             'aquamarine4','darkseagreen','palegreen4',
                             'wheat','tan','dimgray' )) +
  ggtitle('Diet') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(legend.position="none")
dev.off()

pdf(file='figs_diet_and_microbiome_merged/unifrac_dietType2017-2020.pdf', height=5, width=6)
ordinate(micro_rare, "PCoA", "unifrac") %>% 
  plot_ordination(micro_rare, ., color='Diet_type', title = "UniFrac", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Diet type"))+
  scale_color_manual(values=c('#CC79A7', '#E69F00', '#56B4E9')) + geom_point(size=3) +
  theme_bw() + stat_ellipse(level=0.5, lty=2, lwd=0.6) +
  theme(axis.text = element_text(size = 11)) 
dev.off()

ordinate(micro_rare, "PCoA", "bray") %>% 
  plot_ordination(micro_rare, ., color='Diet_type', title = "Bray-Curtis", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Diet type"))+
  scale_color_manual(values=c('#CC79A7', '#E69F00', '#56B4E9')) + geom_point(size=3) +
  theme_bw() + stat_ellipse(level=0.5, lty=2, lwd=0.6) +
  theme(axis.text = element_text(size = 11)) 


#pdf(file='~/Documents/Conferences:Presentations/AOS-SCO2021/hist_dietType.pdf', height=3, width=4)
hist(index_ab$i2, breaks=10, xlab='Index of diet specialization', main='', col='gray87')
dev.off()

md$observed<-estimate_richness(micro_rare, split=T, measures=c('Observed'))
kruskal.test(md$shannon$Shannon,md$Diet_type)
kruskal.test(md$observed$Observed,md$Diet_type)

###for topology analyses of 2017-2020 samples merged 
table(diet_rare@sam_data$Species) #15 species
#convert feature table to biom format
library(biomformat);packageVersion("biomformat") #‘1.16.0’
otu_biom<-as(otu_table(diet_rare),"matrix")
otu_biom<-make_biom(data=otu_biom)
write_biom(otu_biom,"diet_rare_biom.biom")

#write tree to .tre file for qiime import
write.tree(diet_rare@phy_tree, file = "diet_rare.tre", append = FALSE, digits = 10, tree.names = FALSE)

#see commands_qiime_diversity_grouped.txt for qiime stuff, return here to convert distance matrices to UPGMA clustered dendrograms
library(ctc)
#read in grouped distance matrices
by<-read.csv('../grouped_diet/exported_braycurtis/distance-matrix.tsv', sep='\t', header=T, row.names=1)
uy<-read.csv('../grouped_diet/exported_unifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
wy<-read.csv('../grouped_diet/exported_weightedunifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
jy<-read.csv('../grouped_diet/exported_jaccard/distance-matrix.tsv', sep='\t', header=T, row.names=1)

#transform to matrix
by<-as.matrix(by)
by<-as.dist(by, diag=T, upper=F) #dist object for hclust
uy<-as.matrix(uy)
uy<-as.dist(uy, diag=T, upper=F)
wy<-as.matrix(wy)
wy<-as.dist(wy, diag=T, upper=F)
jy<-as.matrix(jy)
jy<-as.dist(jy, diag=T, upper=F)


#create UPGMA dendrograms by clustering
bd<-hclust(by, method='average')
ud<-hclust(uy, method='average')
wd<-hclust(wy, method='average')
jd<-hclust(jy, method='average')
plot(bd, hang=-1)
plot(ud, hang=-1)
plot(wd, hang=-1)
plot(jd, hang=-1)

#write to newick format for TreeCmp
bn<-hc2Newick(bd, flat=T)
un<-hc2Newick(ud, flat=T)
wn<-hc2Newick(wd, flat=T)
jn<-hc2Newick(jd, flat=T)
write.table(bn, '../grouped_diet/topology_analyses/newick_bray_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(un, '../grouped_diet/topology_analyses/newick_unifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(wn, '../grouped_diet/topology_analyses/newick_weightedunifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(jn, '../grouped_diet/topology_analyses/newick_jaccard_grouped.newick', quote=F, col.names=F, row.names=F)


#find subsets for ALL NY
dim(diet_rare@sam_data);dim(micro_rare@sam_data) #same 216 individuals
table(diet_rare@sam_data$State) #135 in NY, 81 in PA
table(diet_rare@sam_data$State, diet_rare@sam_data$Year)
table(diet_rare@sam_data$State, diet_rare@sam_data$Year, diet_rare@sam_data$Species)

barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$State=='NY')]), las=2) #13 spp. with 7+ indiv. each.
ny<-rownames(diet_rare@sam_data)[which(diet_rare@sam_data$State=='NY')] #all 135 NY individuals
write.csv(ny, '../grouped_diet/topology_analyses/subset_NY_merged/individuals_ny.txt', quote=F, row.names=F)
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$State=='NY' & diet_rare@sam_data$Year!='2020')]), las=2) #excludes 2020, still 5+ individuals in NY
ny2<-rownames(diet_rare@sam_data)[which(diet_rare@sam_data$State=='NY'& diet_rare@sam_data$Year!='2020')] #all NY individuals from 2017-2020
write.csv(ny2, '../grouped_diet/topology_analyses/subset_NY_2017to2019/individuals_ny.txt', quote=F, row.names=F)
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$Year=='2020')]), las=2) # 14 spp. with 4+ indiv ea.
tw<-rownames(diet_rare@sam_data)[which(diet_rare@sam_data$Year=='2020')] #all 99 2020 individuals
write.csv(tw, '../grouped_diet/topology_analyses/subset_2020/individuals_2020.txt', quote=F, row.names=F)
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$State=='PA' & diet_rare@sam_data$Year=='2020')]), las=2) #12 spp. More individuals, but singleton spp.
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$State=='NY' & diet_rare@sam_data$Year=='2020')]), las=2) # 12 spp. At least 2 individuals per spp. Use NY 2020.
tn<-rownames(diet_rare@sam_data)[which(diet_rare@sam_data$State=='NY' & diet_rare@sam_data$Year=='2020')] #all 46 NY 2020 individuals
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$State=='NY' & diet_rare@sam_data$Year=='2019')]), las=2) #lots of singleton spp.
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$Year=='2018')]), las=2) #lots of singleton spp, don't use
barplot(table(diet_rare@sam_data$Species[which(diet_rare@sam_data$Year!='2020')]), las=2) #excludes 2020, 6+ indiv/spp.14 spp.
all1<-rownames(diet_rare@sam_data)[which(diet_rare@sam_data$Year!='2020')] #all  individuals from 2017-2019
write.csv(all1, '../grouped_diet/topology_analyses/subset_all1/individuals_all1.txt', quote=F, row.names=F)


#see commands_qiime_diversity_grouped_subset2020.txt for qiime stuff, return here to convert distance matrices to UPGMA clustered dendrograms
library(ctc)
#read in grouped distance matrices
by<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_braycurtis/distance-matrix.tsv', sep='\t', header=T, row.names=1)
uy<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_unifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
wy<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_weightedunifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
jy<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_jaccard/distance-matrix.tsv', sep='\t', header=T, row.names=1)

#transform to matrix
by<-as.matrix(by)
by<-as.dist(by, diag=T, upper=F) #dist object for hclust
uy<-as.matrix(uy)
uy<-as.dist(uy, diag=T, upper=F)
wy<-as.matrix(wy)
wy<-as.dist(wy, diag=T, upper=F)
jy<-as.matrix(jy)
jy<-as.dist(jy, diag=T, upper=F)


#create UPGMA dendrograms by clustering
bd<-hclust(by, method='average')
ud<-hclust(uy, method='average')
wd<-hclust(wy, method='average')
jd<-hclust(jy, method='average')
plot(bd, hang=-1)
plot(ud, hang=-1)
plot(wd, hang=-1)
plot(jd, hang=-1)

#write to newick format for TreeCmp
bn<-hc2Newick(bd, flat=T)
un<-hc2Newick(ud, flat=T)
wn<-hc2Newick(wd, flat=T)
jn<-hc2Newick(jd, flat=T)
write.table(bn, '../grouped_diet/topology_analyses/subset_2020/newick_bray_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(un, '../grouped_diet/topology_analyses/subset_2020/newick_unifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(wn, '../grouped_diet/topology_analyses/subset_2020/newick_weightedunifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(jn, '../grouped_diet/topology_analyses/subset_2020/newick_jaccard_grouped.newick', quote=F, col.names=F, row.names=F)

##subset 2020 microbiome
#read in grouped distance matrices
by<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_braycurtis_micro/distance-matrix.tsv', sep='\t', header=T, row.names=1)
uy<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_unifrac_micro/distance-matrix.tsv', sep='\t', header=T, row.names=1)
wy<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_weightedunifrac_micro/distance-matrix.tsv', sep='\t', header=T, row.names=1)
jy<-read.csv('../grouped_diet/topology_analyses/subset_2020/exported_jaccard_micro/distance-matrix.tsv', sep='\t', header=T, row.names=1)

#transform to matrix
by<-as.matrix(by)
by<-as.dist(by, diag=T, upper=F) #dist object for hclust
uy<-as.matrix(uy)
uy<-as.dist(uy, diag=T, upper=F)
wy<-as.matrix(wy)
wy<-as.dist(wy, diag=T, upper=F)
jy<-as.matrix(jy)
jy<-as.dist(jy, diag=T, upper=F)


#create UPGMA dendrograms by clustering
bd<-hclust(by, method='average')
ud<-hclust(uy, method='average')
wd<-hclust(wy, method='average')
jd<-hclust(jy, method='average')
plot(bd, hang=-1)
plot(ud, hang=-1)
plot(wd, hang=-1)
plot(jd, hang=-1)

#write to newick format for TreeCmp
bn<-hc2Newick(bd, flat=T)
un<-hc2Newick(ud, flat=T)
wn<-hc2Newick(wd, flat=T)
jn<-hc2Newick(jd, flat=T)
write.table(bn, '../grouped_diet/topology_analyses/subset_2020/newick_bray_grouped_micro.newick', quote=F, col.names=F, row.names=F)
write.table(un, '../grouped_diet/topology_analyses/subset_2020/newick_unifrac_grouped_micro.newick', quote=F, col.names=F, row.names=F)
write.table(wn, '../grouped_diet/topology_analyses/subset_2020/newick_weightedunifrac_grouped_micro.newick', quote=F, col.names=F, row.names=F)
write.table(jn, '../grouped_diet/topology_analyses/subset_2020/newick_jaccard_grouped_micro.newick', quote=F, col.names=F, row.names=F)

###subset NY diet
#read in grouped distance matrices
by<-read.csv('../grouped_diet/topology_analyses/subset_NY_2017to2019/exported_braycurtis/distance-matrix.tsv', sep='\t', header=T, row.names=1)
uy<-read.csv('../grouped_diet/topology_analyses/subset_NY_2017to2019/exported_unifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
wy<-read.csv('../grouped_diet/topology_analyses/subset_NY_2017to2019/exported_weightedunifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
jy<-read.csv('../grouped_diet/topology_analyses/subset_NY_2017to2019/exported_jaccard/distance-matrix.tsv', sep='\t', header=T, row.names=1)

#transform to matrix
by<-as.matrix(by)
by<-as.dist(by, diag=T, upper=F) #dist object for hclust
uy<-as.matrix(uy)
uy<-as.dist(uy, diag=T, upper=F)
wy<-as.matrix(wy)
wy<-as.dist(wy, diag=T, upper=F)
jy<-as.matrix(jy)
jy<-as.dist(jy, diag=T, upper=F)

#create UPGMA dendrograms by clustering
bd<-hclust(by, method='average')
ud<-hclust(uy, method='average')
wd<-hclust(wy, method='average')
jd<-hclust(jy, method='average')
plot(bd, hang=-1)
plot(ud, hang=-1)
plot(wd, hang=-1)
plot(jd, hang=-1)

#write to newick format for TreeCmp
bn<-hc2Newick(bd, flat=T)
un<-hc2Newick(ud, flat=T)
wn<-hc2Newick(wd, flat=T)
jn<-hc2Newick(jd, flat=T)
write.table(bn, '../grouped_diet/topology_analyses/subset_NY_2017to2019/newick_bray_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(un, '../grouped_diet/topology_analyses/subset_NY_2017to2019/newick_unifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(wn, '../grouped_diet/topology_analyses/subset_NY_2017to2019/newick_weightedunifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(jn, '../grouped_diet/topology_analyses/subset_NY_2017to2019/newick_jaccard_grouped.newick', quote=F, col.names=F, row.names=F)

###subset all1 (2017-2019) diet
#read in grouped distance matrices
by<-read.csv('../grouped_diet/topology_analyses/subset_all1/exported_braycurtis/distance-matrix.tsv', sep='\t', header=T, row.names=1)
uy<-read.csv('../grouped_diet/topology_analyses/subset_all1/exported_unifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
wy<-read.csv('../grouped_diet/topology_analyses/subset_all1/exported_weightedunifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
jy<-read.csv('../grouped_diet/topology_analyses/subset_all1/exported_jaccard/distance-matrix.tsv', sep='\t', header=T, row.names=1)

#transform to matrix
by<-as.matrix(by)
by<-as.dist(by, diag=T, upper=F) #dist object for hclust
uy<-as.matrix(uy)
uy<-as.dist(uy, diag=T, upper=F)
wy<-as.matrix(wy)
wy<-as.dist(wy, diag=T, upper=F)
jy<-as.matrix(jy)
jy<-as.dist(jy, diag=T, upper=F)

#create UPGMA dendrograms by clustering
bd<-hclust(by, method='average')
ud<-hclust(uy, method='average')
wd<-hclust(wy, method='average')
jd<-hclust(jy, method='average')
plot(bd, hang=-1)
plot(ud, hang=-1)
plot(wd, hang=-1)
plot(jd, hang=-1)

#write to newick format for TreeCmp
bn<-hc2Newick(bd, flat=T)
un<-hc2Newick(ud, flat=T)
wn<-hc2Newick(wd, flat=T)
jn<-hc2Newick(jd, flat=T)
write.table(bn, '../grouped_diet/topology_analyses/subset_all1/newick_bray_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(un, '../grouped_diet/topology_analyses/subset_all1/newick_unifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(wn, '../grouped_diet/topology_analyses/subset_all1/newick_weightedunifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(jn, '../grouped_diet/topology_analyses/subset_all1/newick_jaccard_grouped.newick', quote=F, col.names=F, row.names=F)


#mantel tests
mantel(bray_micro, bray_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.05884, p=0.045
aa = as.vector(bray_micro)
tt = as.vector(bray_diet)
mat = data.frame(aa,tt)

pdf(file='figs_diet_and_microbiome_merged/mantel_dietBC2017-2020.pdf', height=5, width=5)
plot(mat$aa, mat$tt, col=alpha('black', 0.4), pch=16,  ylim=c(0,1), xlim=c(0,1),
     xlab='Microbiome similarity (Bray-Curtis)', ylab='Diet similarity (Bray-Curtis)', cex=1.3)
abline(lm(mat$tt~mat$aa), xpd=F, lty=2, lwd=1.5, col='red')
dev.off()

mantel(unifrac_micro, unifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.1632 , p=0.001
mat$bb = as.vector(unifrac_micro)
mat$cc = as.vector(unifrac_diet)

pdf(file='figs_diet_and_microbiome_merged/mantel_dietUni2017-2020.pdf', height=5, width=5)
plot(mat$bb, mat$cc, col=alpha('black', 0.4), pch=16, ylim=c(0.3,1), xlim=c(0.3,1),
     xlab='Microbiome similarity (UniFrac)', ylab='Diet similarity (UniFrac)', cex=1.3)
abline(lm(mat$cc~mat$bb), xpd=F, lty=2, lwd=1.5, col='red')
dev.off()

wunifrac_micro<-distance(micro_rare,'wunifrac')
wunifrac_diet<-distance(diet_rare,'wunifrac')
mantel(wunifrac_micro, wunifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.02512, p=0.352
mat$dd = as.vector(wunifrac_micro)
mat$ee = as.vector(wunifrac_diet)
plot(mat$dd, mat$ee, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (weighted UniFrac)', ylab='Diet similarity (weighted UniFrac)', cex=1.3)
abline(lm(mat$ee~mat$dd), xpd=F, lty=2, lwd=1.5, col='red')

mantel(jac_micro, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.1744, p=0.001
mat$ff = as.vector(jac_micro)
mat$gg = as.vector(jac_diet)
plot(mat$ff, mat$gg, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (Jaccard)', ylab='Diet similarity (Jaccard)', cex=1.3)
abline(lm(mat$gg~mat$ff), xpd=F, lty=2, lwd=1.5, col='red')

#export distance matricies 
saveRDS(bray_micro, file='bray_micro_full.rds')
saveRDS(jac_micro, file='jac_micro_full.rds')
saveRDS(unifrac_micro, file='unifrac_micro_full.rds')
saveRDS(wunifrac_micro, file='wunifrac_micro_full.rds')
#but don't use these since it is the diet subset--redo using phyloseq.R objects!!
saveRDS(bray_diet, file='bray_diet_full.rds')
saveRDS(jac_diet, file='jac_diet_full.rds')
saveRDS(unifrac_diet, file='unifrac_diet_full.rds')
saveRDS(wunifrac_diet, file='wunifrac_diet_full.rds')

#mantel tests for 2020-subset
#distances for 2020 subset
bray_micro_b2<-distance(subset_samples(micro_rare, Year=='2020'), method="bray")
bray_diet_b2<-distance(subset_samples(diet_rare, Year=='2020'), method="bray")
jac_micro_b2<-distance(subset_samples(micro_rare, Year=='2020'), method="jaccard", binary=T)
jac_diet_b2<-distance(subset_samples(diet_rare, Year=='2020'), method="jaccard", binary=T)
uni_micro_b2<-distance(subset_samples(micro_rare, Year=='2020'), method="unifrac")
uni_diet_b2<-distance(subset_samples(diet_rare, Year=='2020'), method="unifrac")
wuni_micro_b2<-distance(subset_samples(micro_rare, Year=='2020'), method="wunifrac")
wuni_diet_b2<-distance(subset_samples(diet_rare, Year=='2020'), method="wunifrac")

mantel(bray_micro_b2, bray_diet_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.008463, p=0.454
hh = as.vector(bray_micro_b2)
ii = as.vector(bray_diet_b2)
mat2 = data.frame(hh,ii)
plot(mat2$hh, mat2$ii, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (Bray-Curtis)', ylab='Diet similarity (Bray-Curtis)', cex=1.3)
abline(lm(mat2$ii~mat2$hh), xpd=F, lty=2, lwd=1.5, col='red')

mantel(jac_micro_b2, jac_diet_b2, method = "spearman", permutations = 999, na.rm = TRUE) #-0.0293, p=0.687
mat2$jj = as.vector(jac_micro_b2)
mat2$kk = as.vector(jac_diet_b2)
plot(mat2$jj, mat2$kk, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (Jaccard)', ylab='Diet similarity (Jaccard)', cex=1.3)
abline(lm(mat2$kk~mat2$jj), xpd=F, lty=2, lwd=1.5, col='red')

mantel(uni_micro_b2, uni_diet_b2, method = "spearman", permutations = 999, na.rm = TRUE) #-0.01635, p=0.588
mat2$ll = as.vector(uni_micro_b2)
mat2$mm = as.vector(uni_diet_b2)
plot(mat2$ll, mat2$mm, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (UniFrac)', ylab='Diet similarity (UniFrac)', cex=1.3)
abline(lm(mat2$mm~mat2$ll), xpd=F, lty=2, lwd=1.5, col='red')

mantel(wuni_micro_b2, wuni_diet_b2, method = "spearman", permutations = 999, na.rm = TRUE) #0.05623, p=0.237
mat2$nn = as.vector(wuni_micro_b2)
mat2$oo = as.vector(wuni_diet_b2)
plot(mat2$nn, mat2$oo, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (weighted UniFrac)', ylab='Diet similarity (weighted UniFrac)', cex=1.3)
abline(lm(mat2$oo~mat2$nn), xpd=F, lty=2, lwd=1.5, col='red')

#export distance matricies 
saveRDS(bray_micro_b2, file='bray_micro_b2.rds')
saveRDS(jac_micro_b2, file='jac_micro_b2.rds')
saveRDS(uni_micro_b2, file='unifrac_micro_b2.rds')
saveRDS(wuni_micro_b2, file='wunifrac_micro_b2.rds')
#but don't use these since it is the diet subset--redo using phyloseq.R objects!!
saveRDS(bray_diet_b2, file='bray_diet_b2.rds')
saveRDS(jac_diet_b2, file='jac_diet_b2.rds')
saveRDS(uni_diet_b2, file='unifrac_diet_b2.rds')
saveRDS(wuni_diet_b2, file='wunifrac_diet_b2.rds')

###diet index, correlatin between 2 data sets
pal<-palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
               "red2", "#75DED6", "#75E563", "black", "#B44AE1"))
index_ab_batch1<-read.csv('~/Documents/REU2021/warbler_data/index_ab.csv')
pdf(file='figs_diet_and_microbiome_merged/index_correlation.pdf', height=5, width=5)
par(mar=c(5, 4, 2, 6.5), xpd=TRUE)
plot(index_ab$i1[1:14], index_ab_batch1$i1, ylim=c(2,2.9), xlim=c(2,2.9),
     col=factor(index_ab_batch1$X), pch=16, cex=1.3,
     ylab='Diet index--batch 1', xlab='Diet index--full dataset') #leave off WEWA
abline(lm(index_ab_batch1$i1~index_ab$i1[1:14]), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(index_ab_batch1$X)), col=1:nlevels(factor(index_ab_batch1$X)), pch=16, inset=c(-0.37,0), pt.cex=1.2,bty='n')
dev.off()
cor.test(index_ab$i1[1:14], index_ab_batch1$i1, method='kendall') #tau=0.5164835, p=0.009753

#also test excluding 4 spp that differ btwn batch 1 and full
include<-c('MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','CAWA','NAWA','OVEN')
kruskal.test(md$shannon$Shannon[which(md$Species %in% include)],md$Diet_type[which(md$Species %in% include)])

#diet index with other new alpha and beta metrics
#original calculate indexes
#i1 is shannon+BC
#i2 is shannon+uni
#calculate chao1 and faith pd for diet alpha
alpha_diet$chao<-estimate_richness(diet_rare, measures=c('Chao1'))[1] #add to alpha dataframe
library(btools) #for faiths pd
testpd<-estimate_pd(diet_rare)
alpha_diet$PD<-testpd$PD
#correlation between alphas
plot(alpha_diet$Shannon~alpha_diet$chao$Chao1)
plot(alpha_diet$Shannon~alpha_diet$PD)
plot(alpha_diet$chao$Chao1~alpha_diet$PD)
#for supplemental table
cor.test(alpha_diet$Shannon, alpha_diet$PD) #pearson r=0.5197285, p<0.001
cor.test(alpha_diet$Shannon, alpha_diet$chao$Chao1) #r=0.4779385, p<0.001
cor.test(alpha_diet$PD, alpha_diet$chao$Chao1) #r=0.7903934, p<0.001

#Jaccard distance
jac_diet<-distance(diet_rare, method='jaccard', binary=T)
jac_micro<-distance(micro_rare, method='jaccard', binary=T)

wunifrac_diet<-distance(diet_rare,'wunifrac')
wunifrac_micro<-distance(micro_rare,'wunifrac')

#for supplmental table
mantel(bray_diet, unifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.3076, p=0.001
mantel(bray_diet, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.6221, p=0.001
mantel(unifrac_diet, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.6346, p=0.001
mantel(bray_diet, wunifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.1383, p=0.001
mantel(jac_diet, wunifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.09229, p=0.014
mantel(unifrac_diet, wunifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.04086, p=0.186

##correlation between alpha diet and alpha micro w/pd and chao
div$diet_PD<-testpd$PD
div$diet_chao1<-alpha_diet$chao$Chao1
#calculate micro pd and chao
alpha_micro$chao<-estimate_richness(micro_rare, measures=c('Chao1'))[1] #add to alpha dataframe
testpd2<-estimate_pd(micro_rare)
alpha_micro$PD<-testpd2$PD
div$micro_PD<-testpd2$PD
div$micro_chao1<-alpha_micro$chao$Chao1
plot(div$micro_PD, div$diet_PD, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (PD)', ylab='Diet alpha diversity (PD)',
     main='2017-2020')
abline(lm(div$diet_PD~div$micro_PD), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(div$micro_PD, div$diet_PD, method='kendall') #tau=-0.03385013, p=0.4592
plot(div$micro_chao1, div$diet_chao1, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (Chao1)', ylab='Diet alpha diversity (Chao1)',
     main='2017-2020')
abline(lm(div$diet_chao1~div$micro_chao1), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(div$micro_chao1, div$diet_chao1, method='kendall') #tau=-0.06261447, p=0.1712

#adonis on diet type
adonis2(jac_micro~Diet_type, data=md) #r2=0.00998, p=0.167

#individual relative abundance plot
level_order<-c('BTNW','MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','HOWA','CAWA','COYE','NAWA','BAWW','WEWA','OVEN')
level_order2<-data.frame(id=mps1$Sample,species=mps1$sample_Species)
orders<-data.frame(species=level_order, order=1:15)
#merge <- left_join(orders, level_order2, by = "species")
#merge$Sample<-merge$id
#orders[14,2]<-14
mps1$order<-match(mps1$sample_Species, orders$species)
mps1$order2<-paste(mps1$order, mps1$Sample, sep='_')
mps1$new<-mps1$order
mps1$new <- car::recode(mps1$new, '10=91;11=92;12=93;13=94;14=95;15=96')
table(mps1$new) 
table(mps1$sample_Species)
mps1$order3<-paste(mps1$new, mps1$Sample, sep='_')

pdf(file='~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/figs/ind_abund_diet.pdf', height=4, width=12)
ggplot(mps1, aes(x = order3, y = Abundance, fill = Order2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','aquamarine3',
                             'aquamarine4','darkseagreen','palegreen4',
                             'wheat','tan','dimgray' )) +
  ggtitle('COI') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(axis.text = element_text(size = 11)) +
  theme(legend.position="none") + geom_vline(xintercept=23.5 ) + geom_vline(xintercept=37.5 ) +
  geom_vline(xintercept=52.5 ) + geom_vline(xintercept=69.5 ) + geom_vline(xintercept= 81.5) + 
  geom_vline(xintercept=93.5 ) + geom_vline(xintercept=113.5 ) + geom_vline(xintercept=129.5 ) +
  geom_vline(xintercept=140.5 ) + geom_vline(xintercept=155.5 ) + geom_vline(xintercept=168.5 ) +
  geom_vline(xintercept=175.5 ) + geom_vline(xintercept=193.5 ) + geom_vline(xintercept=198.5 ) +
  geom_vline(xintercept=216.5) + theme(axis.text.x = element_blank())
dev.off()

table(md$Species)

#kruskal test of gm by diet type for chao and pd
kruskal.test(alpha_micro$chao$Chao1, sample_data(micro_rare)$Diet_type) #ns chi-squared = 0.31338, df = 2, p-value = 0.855
kruskal.test(alpha_micro$PD, sample_data(micro_rare)$Diet_type) #ns chi-squared = 0.1334, df = 2, p-value = 0.9355

#extremes<-as.data.frame(div), from above
#extremes$Diet_type<-md$Diet_type
#extremes<-extremes[which(!(extremes$Diet_type=='Intermediate')),] #93 individuals, just generalists and specailists
kruskal.test(extremes$micro_chao1, extremes$Diet_type) #ns chi-squared = 0.28767, df = 1, p-value = 0.5917
kruskal.test(extremes$micro_PD, extremes$Diet_type) #ns chi-squared = 0.17475, df = 1, p-value = 0.6759

#PGLS
wtre<-read.tree("~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/eliot_example/scaledWarbler_full.tre")
wtre$tip.label<-c('OVEN','WEWA','BAWW','NAWA','COYE','HOWA','AMRE','BTBW','NOPA','MAWA','BLBW','CSWA','MYWA','BTNW','CAWA')
wtre$root.edge<-0
library(phylolm)
#cor.test was (div$micro_chao1, div$diet_chao1)

fit_shan = phylolm(alpha_micro~alpha_diet,data=index_ab,phy=wtre,model="lambda", boot=1000) #micro is dependent variable (does micro depend on diet?)
summary(fit_shan) #coef=-0.40536, p=0.2729
#paper that used median species values in phylolm https://onlinelibrary-wiley-com.ezaccess.libraries.psu.edu/doi/10.1111/ibi.12420

#add species averages for other alpha metrics
index_ab$diet_chao1=c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15)
index_ab$micro_chao1=c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
index_ab$diet_pd=c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15)
index_ab$micro_pd=c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15)

fit_chao = phylolm(micro_chao1~diet_chao1,data=index_ab,phy=wtre,model="lambda", boot=1000) #micro is dependent variable (does micro depend on diet?)
summary(fit_chao) #coef=-0.023670, p=0.955196

fit_pd = phylolm(micro_pd~diet_pd,data=index_ab,phy=wtre,model="lambda", boot=1000) #micro is dependent variable (does micro depend on diet?)
summary(fit_pd) #coef=0.33992, p=0.62388

#standard deviation for N individuals per species in these analyses
sd(table(micro_rare@sam_data$Species)) #4.7
