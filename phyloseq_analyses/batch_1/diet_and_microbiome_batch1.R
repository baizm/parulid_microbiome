#script for batch 1 diet & microbiome analyses with the tree output from amptk, rooted here on random arthropod
setwd("~/Documents/REU2021/warbler_data/")
load('diet_and_microbiome2.RData')
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
library(randomcoloR)

#read in micro_rare2 from old script, remove zero-sum taxa. Here there are no spp with < 5 individuals.
micro_rare2<-readRDS('micro_rare2.rds')
micro_rare2<-prune_taxa(taxa_sums(micro_rare2) > 0, micro_rare2) #130 samples, 4471 taxa

#read in diet data, root tree, do diet & microbiome analyses
####------read in OTU table, has to be a matrix---------------------
otu_table_diet<-read.csv('co1_feature-table.tsv', sep='\t', header=T)
otu_table_diet<-otu_table_diet[,order(colnames(otu_table_diet))]
#convert to matrix
otumat_diet<-as.matrix(otu_table_diet)
###-------read in taxonomy table-----------------
tax_table_diet<-read.csv('co1_taxonomy.tsv', sep='\t', header=T)
dim(tax_table_diet) #3949 OTUs total before filtering
#filter out non-arthropoda
tax_table2_diet<-tax_table_diet[which(tax_table_diet$Phylum=='p:Arthropoda'),] #note: analagous to 2017-2020 merged since there are no unassigned OTUs here, just excludes chordates. Run table(tax_table_diet$Phylum, useNA = 'always')
#convert to matrix
taxmat_diet<-as.matrix(tax_table2_diet)
#filter out non-arthropoda from OTU table
otumat_diet<-otumat_diet[which(rownames(otumat_diet) %in% rownames(taxmat_diet)),]
###-------read in sample data---------------------
sampledata<-read.csv('~/Documents/Toews_Lab/16s_data/Jan2020/metadata2_16s_ny.tsv', sep='\t', header=T)
rownames(sampledata)<-sampledata$id
x<-sampledata$id[(which(!(sampledata$id %in% colnames(otumat_diet))))]
sampledata<-sampledata[which(!(sampledata$id %in% x)),] #get rid of 3 samples filtered by qiime (2 neg + 1 low)
sampledata<-sampledata[,2:29]
###-------read in nwk tree with ape package-----------
#tree_diet<-read.tree('~/Documents/Toews_Lab/warbler_diet/2017-2019_Marcella/out.tree.phy') #used the one below, not sure if they differ, but andrew ran the first data
tree_diet<-read.tree('all_samples.tree.phy') #this is from Andrew, instead of co1_trees/clustered.filtered.otus.curated.50plusreads.NJ.nh, like the old script
###----------combine into phyloseq object-------------
OTU_diet = otu_table(otumat_diet, taxa_are_rows = TRUE)
TAX_diet = tax_table(taxmat_diet)
SAM = sample_data(sampledata)
physeq_diet = phyloseq(OTU_diet, TAX_diet, SAM, tree_diet)

sample_sums(physeq_diet)

#rarefy diet object and filter to individuals in micro_rare2
hist(colSums(physeq_diet@otu_table))
summary(colSums(physeq_diet@otu_table))
rarecurve(t(otu_table(physeq_diet)), step=50, lwd=0.8, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,15000), ylim=c(0,200))

diet_rare2<-prune_samples(rownames(physeq_diet@sam_data)[which(rownames(physeq_diet@sam_data) %in% rownames(micro_rare2@sam_data))], physeq_diet)
summary(sample_sums(diet_rare2))
diet_rare2<-rarefy_even_depth(diet_rare2, sample.size=8500, rngseed=8, replace=F)


rownames(micro_rare2@sam_data)[which(!(rownames(micro_rare2@sam_data) %in% rownames(diet_rare2@sam_data)))] #all there!

#root the diet tree!
arachnids<-rownames(diet_rare2@tax_table)[which(diet_rare2@tax_table[,3]=='c:Arachnida')] #build vector of arachnid OTUs
sample(arachnids, 1) #choose random arachnid to root the tree on, got OTU621 (OTU1598 not present, was used in 2017-2020 dataset, but sure if names are consistent though)
diet_rare2@phy_tree<-root(diet_rare2@phy_tree, outgroup = 'OTU621', resolve.root=T)

table(diet_rare2@sam_data$Species) #none fewer than 6 indiv.
diet_rare2<-prune_taxa(taxa_sums(diet_rare2) > 0, diet_rare2) #130 samples, 1773 taxa

ordinate(diet_rare2, "PCoA", "bray") %>% 
  plot_ordination(diet_rare2, ., color='Species', title = "Bray-Curtis", type='samples') +
  geom_point(size=2)

ps0<-tax_glom(diet_rare2, taxrank = "Class", NArm=F) #glom by class, not helpful
mps0<-psmelt(ps0)
table(mps0$Class)
#pdf(file='taxaPlot_rarefied.pdf', height=6, width=8)
ggplot(mps0, aes(x = sample_Species, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('cornsilk3','darkblue',
                             'cornflowerblue','coral3','darkgreen',
                             'darkolivegreen3','black','red')) +
  ggtitle('Diet') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
#dev.off()

###diet relative abundance plot distinguishing only top 7 orders
topc<-as.data.frame(diet_rare2@tax_table)
sort(table(topc$Order), decreasing=T)
topcc<-names(sort(table(topc$Order), decreasing=T)[1:7]) #top 7 orders

ps1<-tax_glom(diet_rare2, taxrank = "Order", NArm=F) #glom by order
mps1<-psmelt(ps1)
mps1$Order[which(!(mps1$Order %in% topcc))]<-'other' #change non-common orders to other 
table(mps1$Order)

mps1$Order2<-as.factor(mps1$Order)
mps1$Order2<-relevel(mps1$Order2, 'other')

mps1$Order[which(mps1$Class=='c:Arachnida' & mps1$Abundance>0)]

ggplot(mps1, aes(x = sample_Species, y = Abundance, fill = Order2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','aquamarine3',
                             'aquamarine4','darkseagreen','palegreen4',
                             'wheat','tan','dimgray' )) +
  ggtitle('Diet') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) 
 # theme(legend.position="none")


#6: calculate alpha and beta for each dataset
alpha_micro<-diversity(micro_rare2@otu_table, index='shannon', MARGIN=2)
alpha_diet<-diversity(diet_rare2@otu_table, index='shannon', MARGIN=2)
plot(alpha_diet, alpha_micro)
abline(lm(alpha_micro~alpha_diet))

div<-micro_rare2@sam_data
div$micro_shannon<-alpha_micro
div$diet_shannon<-alpha_diet

library(randomcoloR)
#palette(distinctColorPalette(14))
palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
          "red2", "#75DED6", "#75E563", "black", "#B44AE1"))

pdf(file='figs_diet_and_microbiome2/diet_micro_scatter2017-2019.pdf', height=5, width=5)
par(mar=c(5, 4, 2, 7.3), xpd=TRUE)
pal<-palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
               "red2", "#75DED6", "#75E563", "black", "#B44AE1"))
plot(div$micro_shannon, div$diet_shannon, col=factor(div$Species), pch=16,
     xlab='Microbiome alpha diversity (Shannon)', ylab='Diet alpha diversity (Shannon)', cex=1.3)
abline(lm(div$diet_shannon~div$micro_shannon), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=16, inset=c(-0.37,0), pt.cex=1.2,bty='n')
dev.off()

alpha_micro<-estimate_richness(micro_rare2, split=T, measures=c('Observed', 'Shannon'))
alpha_diet<-estimate_richness(diet_rare2, split=T, measures=c('Observed', 'Shannon'))

plot(alpha_diet$Observed, alpha_micro$Observed)
abline(lm(alpha_micro$Observed~alpha_diet$Observed))

div$observed_micro<-alpha_micro$Observed
div$observed_diet<-alpha_diet$Observed

plot_richness(diet_rare2, x="Species", color="Species",measures="Observed") + theme(legend.position="none") +
  geom_boxplot() 

###--find mean alpha diversity for each species, diet and micro
a1<-mean(div$micro_shannon[which(div$Species=='AMRE')])
b1<-mean(div$diet_shannon[which(div$Species=='AMRE')])
c1<-mean(div$micro_chao1[which(div$Species=='AMRE')])
d1<-mean(div$diet_chao1[which(div$Species=='AMRE')])
e1<-mean(div$micro_PD[which(div$Species=='AMRE')])
f1<-mean(div$diet_PD[which(div$Species=='AMRE')])

a2<-mean(div$micro_shannon[which(div$Species=='BAWW')])
b2<-mean(div$diet_shannon[which(div$Species=='BAWW')])
c2<-mean(div$micro_chao1[which(div$Species=='BAWW')])
d2<-mean(div$diet_chao1[which(div$Species=='BAWW')])
e2<-mean(div$micro_PD[which(div$Species=='BAWW')])
f2<-mean(div$diet_PD[which(div$Species=='BAWW')])

a3<-mean(div$micro_shannon[which(div$Species=='BLBW')])
b3<-mean(div$diet_shannon[which(div$Species=='BLBW')])
c3<-mean(div$micro_chao1[which(div$Species=='BLBW')])
d3<-mean(div$diet_chao1[which(div$Species=='BLBW')])
e3<-mean(div$micro_PD[which(div$Species=='BLBW')])
f3<-mean(div$diet_PD[which(div$Species=='BLBW')])

a4<-mean(div$micro_shannon[which(div$Species=='BTBW')])
b4<-mean(div$diet_shannon[which(div$Species=='BTBW')])
c4<-mean(div$micro_chao1[which(div$Species=='BTBW')])
d4<-mean(div$diet_chao1[which(div$Species=='BTBW')])
e4<-mean(div$micro_PD[which(div$Species=='BTBW')])
f4<-mean(div$diet_PD[which(div$Species=='BTBW')])

a5<-mean(div$micro_shannon[which(div$Species=='BTNW')])
b5<-mean(div$diet_shannon[which(div$Species=='BTNW')])
c5<-mean(div$micro_chao1[which(div$Species=='BTNW')])
d5<-mean(div$diet_chao1[which(div$Species=='BTNW')])
e5<-mean(div$micro_PD[which(div$Species=='BTNW')])
f5<-mean(div$diet_PD[which(div$Species=='BTNW')])

a6<-mean(div$micro_shannon[which(div$Species=='CAWA')])
b6<-mean(div$diet_shannon[which(div$Species=='CAWA')])
c6<-mean(div$micro_chao1[which(div$Species=='CAWA')])
d6<-mean(div$diet_chao1[which(div$Species=='CAWA')])
e6<-mean(div$micro_PD[which(div$Species=='CAWA')])
f6<-mean(div$diet_PD[which(div$Species=='CAWA')])

a7<-mean(div$micro_shannon[which(div$Species=='COYE')])
b7<-mean(div$diet_shannon[which(div$Species=='COYE')])
c7<-mean(div$micro_chao1[which(div$Species=='COYE')])
d7<-mean(div$diet_chao1[which(div$Species=='COYE')])
e7<-mean(div$micro_PD[which(div$Species=='COYE')])
f7<-mean(div$diet_PD[which(div$Species=='COYE')])

a8<-mean(div$micro_shannon[which(div$Species=='CSWA')])
b8<-mean(div$diet_shannon[which(div$Species=='CSWA')])
c8<-mean(div$micro_chao1[which(div$Species=='CSWA')])
d8<-mean(div$diet_chao1[which(div$Species=='CSWA')])
e8<-mean(div$micro_PD[which(div$Species=='CSWA')])
f8<-mean(div$diet_PD[which(div$Species=='CSWA')])

a9<-mean(div$micro_shannon[which(div$Species=='HOWA')])
b9<-mean(div$diet_shannon[which(div$Species=='HOWA')])
c9<-mean(div$micro_chao1[which(div$Species=='HOWA')])
d9<-mean(div$diet_chao1[which(div$Species=='HOWA')])
e9<-mean(div$micro_PD[which(div$Species=='HOWA')])
f9<-mean(div$diet_PD[which(div$Species=='HOWA')])

a10<-mean(div$micro_shannon[which(div$Species=='MAWA')])
b10<-mean(div$diet_shannon[which(div$Species=='MAWA')])
c10<-mean(div$micro_chao1[which(div$Species=='MAWA')])
d10<-mean(div$diet_chao1[which(div$Species=='MAWA')])
e10<-mean(div$micro_PD[which(div$Species=='MAWA')])
f10<-mean(div$diet_PD[which(div$Species=='MAWA')])

a11<-mean(div$micro_shannon[which(div$Species=='MYWA')])
b11<-mean(div$diet_shannon[which(div$Species=='MYWA')])
c11<-mean(div$micro_chao1[which(div$Species=='MYWA')])
d11<-mean(div$diet_chao1[which(div$Species=='MYWA')])
e11<-mean(div$micro_PD[which(div$Species=='MYWA')])
f11<-mean(div$diet_PD[which(div$Species=='MYWA')])

a12<-mean(div$micro_shannon[which(div$Species=='NAWA')])
b12<-mean(div$diet_shannon[which(div$Species=='NAWA')])
c12<-mean(div$micro_chao1[which(div$Species=='NAWA')])
d12<-mean(div$diet_chao1[which(div$Species=='NAWA')])
e12<-mean(div$micro_PD[which(div$Species=='NAWA')])
f12<-mean(div$diet_PD[which(div$Species=='NAWA')])

a13<-mean(div$micro_shannon[which(div$Species=='NOPA')])
b13<-mean(div$diet_shannon[which(div$Species=='NOPA')])
c13<-mean(div$micro_chao1[which(div$Species=='NOPA')])
d13<-mean(div$diet_chao1[which(div$Species=='NOPA')])
e13<-mean(div$micro_PD[which(div$Species=='NOPA')])
f13<-mean(div$diet_PD[which(div$Species=='NOPA')])

a14<-mean(div$micro_shannon[which(div$Species=='OVEN')])
b14<-mean(div$diet_shannon[which(div$Species=='OVEN')])
c14<-mean(div$micro_chao1[which(div$Species=='OVEN')])
d14<-mean(div$diet_chao1[which(div$Species=='OVEN')])
e14<-mean(div$micro_PD[which(div$Species=='OVEN')])
f14<-mean(div$diet_PD[which(div$Species=='OVEN')])

#calculate beta diversity!
bray1a<-mean(distance(subset_samples(micro_rare2, Species=='AMRE'), "bray"))
bray1b<-mean(distance(subset_samples(diet_rare2, Species=='AMRE'), "bray"))
bray2a<-mean(distance(subset_samples(micro_rare2, Species=='BAWW'), "bray"))
bray2b<-mean(distance(subset_samples(diet_rare2, Species=='BAWW'), "bray"))
bray3a<-mean(distance(subset_samples(micro_rare2, Species=='BLBW'), "bray"))
bray3b<-mean(distance(subset_samples(diet_rare2, Species=='BLBW'), "bray"))
bray4a<-mean(distance(subset_samples(micro_rare2, Species=='BTBW'), "bray"))
bray4b<-mean(distance(subset_samples(diet_rare2, Species=='BTBW'), "bray"))
bray5a<-mean(distance(subset_samples(micro_rare2, Species=='BTNW'), "bray"))
bray5b<-mean(distance(subset_samples(diet_rare2, Species=='BTNW'), "bray"))
bray6a<-mean(distance(subset_samples(micro_rare2, Species=='CAWA'), "bray"))
bray6b<-mean(distance(subset_samples(diet_rare2, Species=='CAWA'), "bray"))
bray7a<-mean(distance(subset_samples(micro_rare2, Species=='COYE'), "bray"))
bray7b<-mean(distance(subset_samples(diet_rare2, Species=='COYE'), "bray"))
bray8a<-mean(distance(subset_samples(micro_rare2, Species=='CSWA'), "bray"))
bray8b<-mean(distance(subset_samples(diet_rare2, Species=='CSWA'), "bray"))
bray9a<-mean(distance(subset_samples(micro_rare2, Species=='HOWA'), "bray"))
bray9b<-mean(distance(subset_samples(diet_rare2, Species=='HOWA'), "bray"))
bray10a<-mean(distance(subset_samples(micro_rare2, Species=='MAWA'), "bray"))
bray10b<-mean(distance(subset_samples(diet_rare2, Species=='MAWA'), "bray"))
bray11a<-mean(distance(subset_samples(micro_rare2, Species=='MYWA'), "bray"))
bray11b<-mean(distance(subset_samples(diet_rare2, Species=='MYWA'), "bray"))
bray12a<-mean(distance(subset_samples(micro_rare2, Species=='NAWA'), "bray"))
bray12b<-mean(distance(subset_samples(diet_rare2, Species=='NAWA'), "bray"))
bray13a<-mean(distance(subset_samples(micro_rare2, Species=='NOPA'), "bray"))
bray13b<-mean(distance(subset_samples(diet_rare2, Species=='NOPA'), "bray"))
bray14a<-mean(distance(subset_samples(micro_rare2, Species=='OVEN'), "bray"))
bray14b<-mean(distance(subset_samples(diet_rare2, Species=='OVEN'), "bray"))

#unifrac
unifrac1a<-mean(distance(subset_samples(micro_rare2,Species=='AMRE'),'unifrac'))
unifrac1b<-mean(distance(subset_samples(diet_rare2,Species=='AMRE'),'unifrac'))
unifrac2a<-mean(distance(subset_samples(micro_rare2,Species=='BAWW'),'unifrac'))
unifrac2b<-mean(distance(subset_samples(diet_rare2,Species=='BAWW'),'unifrac'))
unifrac3a<-mean(distance(subset_samples(micro_rare2,Species=='BLBW'),'unifrac'))
unifrac3b<-mean(distance(subset_samples(diet_rare2,Species=='BLBW'),'unifrac'))
unifrac4a<-mean(distance(subset_samples(micro_rare2,Species=='BTBW'),'unifrac'))
unifrac4b<-mean(distance(subset_samples(diet_rare2,Species=='BTBW'),'unifrac'))
unifrac5a<-mean(distance(subset_samples(micro_rare2,Species=='BTNW'),'unifrac'))
unifrac5b<-mean(distance(subset_samples(diet_rare2,Species=='BTNW'),'unifrac'))
unifrac6a<-mean(distance(subset_samples(micro_rare2,Species=='CAWA'),'unifrac'))
unifrac6b<-mean(distance(subset_samples(diet_rare2,Species=='CAWA'),'unifrac'))
unifrac7a<-mean(distance(subset_samples(micro_rare2,Species=='COYE'),'unifrac'))
unifrac7b<-mean(distance(subset_samples(diet_rare2,Species=='COYE'),'unifrac'))
unifrac8a<-mean(distance(subset_samples(micro_rare2,Species=='CSWA'),'unifrac'))
unifrac8b<-mean(distance(subset_samples(diet_rare2,Species=='CSWA'),'unifrac'))
unifrac9a<-mean(distance(subset_samples(micro_rare2,Species=='HOWA'),'unifrac'))
unifrac9b<-mean(distance(subset_samples(diet_rare2,Species=='HOWA'),'unifrac'))
unifrac10a<-mean(distance(subset_samples(micro_rare2,Species=='MAWA'),'unifrac'))
unifrac10b<-mean(distance(subset_samples(diet_rare2,Species=='MAWA'),'unifrac'))
unifrac11a<-mean(distance(subset_samples(micro_rare2,Species=='MYWA'),'unifrac'))
unifrac11b<-mean(distance(subset_samples(diet_rare2,Species=='MYWA'),'unifrac'))
unifrac12a<-mean(distance(subset_samples(micro_rare2,Species=='NAWA'),'unifrac'))
unifrac12b<-mean(distance(subset_samples(diet_rare2,Species=='NAWA'),'unifrac'))
unifrac13a<-mean(distance(subset_samples(micro_rare2,Species=='NOPA'), 'unifrac'))
unifrac13b<-mean(distance(subset_samples(diet_rare2,Species=='NOPA'), 'unifrac'))
unifrac14a<-mean(distance(subset_samples(micro_rare2,Species=='OVEN'),'unifrac'))
unifrac14b<-mean(distance(subset_samples(diet_rare2,Species=='OVEN'),'unifrac'))

#build dataframe for diet index score!
index_ab<-data.frame("alpha_diet"=c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14),
                     "alpha_micro"=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14),
                     "beta_diet"=c(bray1b,bray2b,bray3b,bray4b,bray5b,bray6b,bray7b,bray8b,bray9b,bray10b,bray11b,bray12b,bray13b,bray14b),
                     "beta_micro"=c(bray1a,bray2a,bray3a,bray4a,bray5a,bray6a,bray7a,bray8a,bray9a,bray10a,bray11a,bray12a,bray13a,bray14a),
                     'unifrac_diet'=c(unifrac1b,unifrac2b,unifrac3b,unifrac4b,unifrac5b,unifrac6b,unifrac7b,unifrac8b,unifrac9b,unifrac10b,unifrac11b,unifrac12b,unifrac13b,unifrac14b),
                     'unifrac_micro'=c(unifrac1a,unifrac2a,unifrac3a,unifrac4a,unifrac5a,unifrac6a,unifrac7a,unifrac8a,unifrac9a,unifrac10a,unifrac11a,unifrac12a,unifrac13a,unifrac14a))
rownames(index_ab)<-c('AMRE','BAWW','BLBW','BTBW','BTNW','CAWA','COYE','CSWA','HOWA','MAWA','MYWA','NAWA','NOPA','OVEN')

plot(index_ab$alpha_diet, index_ab$beta_diet)
plot(index_ab$alpha_diet, index_ab$unifrac_diet)
plot(index_ab$alpha_diet, index_ab$alpha_micro)
abline(lm(index_ab$alpha_micro~index_ab$alpha_diet))
summary(lm(index_ab$alpha_micro~index_ab$alpha_diet)) #ns

plot(index_ab$beta_diet, index_ab$unifrac_diet)
max(index_ab$beta_diet)-min(index_ab$beta_diet) #larger span in bray 
max(index_ab$unifrac_diet)-min(index_ab$unifrac_diet)

#calculate indexes
index_ab$i1<-index_ab$alpha_diet+index_ab$beta_diet
index_ab$i2<-index_ab$alpha_diet+index_ab$unifrac_diet
write.csv(index_ab, 'index_ab.csv')

plot(index_ab$i1, index_ab$i2)
plot(index_ab$i2, index_ab$alpha_micro)
abline(lm(index_ab$alpha_micro~index_ab$i2))


#PERMANOVAs
#generate new distance matrices
bray_micro2<-distance(micro_rare2, 'bray')
bray_diet2<-distance(diet_rare2, 'bray')
unifrac_micro<-distance(micro_rare2,'unifrac')
unifrac_diet<-distance(diet_rare2,'unifrac')
weightedu_micro<-distance(micro_rare2, 'wunifrac')
weightedu_diet<-distance(diet_rare2, 'wunifrac')

md<-data.frame(diet_rare2@sam_data, row.names=rownames(diet_rare2@sam_data)) #build from this as micro_rare2@sam_data has diet type from old script
md$Year<-as.factor(md$Year) #make non-continuous

adonis(bray_micro2 ~ Species, data=md) #r2=0.11057, p=0.062
adonis(bray_micro2 ~ Year, data=md) #r2=0.01958, p=0.065
adonis(bray_micro2 ~ State, data=md) #r2=0.02016, p=0.001***
adonis(bray_micro2 ~ Diet_type, data=md) #r2=0.01963,  p=0.062

adonis(unifrac_micro~Species,data=md) #r2=0.12751, p=0.001***
adonis(unifrac_micro~Year,data=md) #r2=0.02286, p=0.002***
adonis(unifrac_micro~State,data=md) #r2=0.01841, p=0.002***
adonis(unifrac_micro~Diet_type, data=md) #r2=0.02189, p=0.013***
permutest(betadisper(unifrac_micro, md$Diet_type)) #ns F=2.1606, p=0.118

adonis(weightedu_micro~Species,data=md) #r2=0.10842, p=0.374
adonis(weightedu_micro~Year,data=md) #r2=0.0226, p=0.21
adonis(weightedu_micro~State,data=md) #r2=0.01827, p=0.082
adonis(weightedu_micro~Diet_type, data=md) #r2=0.02114, p=0.246

#pdf(file='figs_diet_and_microbiome2/hist_index_scores2017-2019.pdf', height=6, width=5)
par(mfrow=c(2,1), mar=c(4,4,2,3))
hist(index_ab$i1, breaks=10, xlab='Index score', main='Index of diet specialization using BC distance')
hist(index_ab$i2, breaks=10, xlab='Index score', main='Index of diet specialization using UniFrac distance')
dev.off()


h <- hist(index_ab$i1, breaks=10, xlab='Index score', plot=F)
cuts <- cut(h$breaks, c(-Inf,2.1,2.6,Inf))
pdf(file='figs_diet_and_microbiome2/hist_colored_scores2017-2019.pdf', height=3.5, width=5)
plot(h, col=alpha(c('gray85', 'gray50', 'gray1')[cuts], 0.7), xlab='Index of diet specialization', main='',cex.lab=1.4)
dev.off()
#'#56B4E9', '#E69F00','#CC79A7'
#
pdf(file='figs_diet_and_microbiome2/plot_index_scores2017-2019.pdf', height=5, width=4)
plot(index_ab$i1~index_ab$i2, ylab='Diet index using Bray-Curtis distance',
     xlab='Diet index using UniFrac distance', main='2017-2019')
dev.off()

index_ab[order(index_ab$i1),]
index_ab[order(index_ab$i2),] #same 2 spp. as generalists and specialists

index_ab$Diet_Type<-c('Specialist','Intermediate','Intermediate',
                      'Intermediate','Generalist','Generalist',
                      'Intermediate','Specialist','Intermediate',
                      'Intermediate','Intermediate','Intermediate',
                      'Intermediate','Intermediate')

md$Diet_type<-index_ab$Diet_Type[match(md$Species, rownames(index_ab))]

micro_rare2@sam_data$Diet_type<-md$Diet_type #corrected diet_type column in micro_rare2!

ordinate(micro_rare2, "PCoA", "bray") %>% 
  plot_ordination(micro_rare2, ., color='Diet_type', title = "Bray-Curtis", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "New Legend Title"))
ordinate(micro_rare2, "PCoA", "unifrac") %>% 
  plot_ordination(micro_rare2, ., color='Diet_type', title = "UniFrac", type='samples') +
  geom_point(size=2)
pdf(file='figs_diet_and_microbiome2/boxplot_shannon_diettype2017-2019.pdf', height=3.5, width=3.5)
plot_richness(micro_rare2, x="Diet_type", color="Diet_type",measures="Shannon") + theme(legend.position="none") +
  geom_boxplot(alpha=0.7, width=0.5) +scale_color_manual(values=c('gray20', 'gray50', 'gray70')) +
  theme_bw() + scale_x_discrete(limits = rev(levels(factor(micro_rare2@sam_data$Diet_type))),labels=c("Low\ndiversity","Intermediate","High\ndiversity")) +
  xlab('Diet type') + ylab('Microbiome alpha diversity (Shannon)') +
  theme(axis.text = element_text(size = 11)) + ylim(0,6.5) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(legend.position="none")
#+scale_color_manual(values=c('#CC79A7', '#E69F00', '#56B4E9'))
dev.off()
plot_richness(micro_rare2, x="Diet_type", color="Diet_type",measures="Observed") + theme(legend.position="none") +
  geom_boxplot(alpha=0.7, width=0.5) +scale_color_manual(values=c('chocolate4', 'aquamarine4', 'darkgoldenrod')) +
  theme_bw()

kruskal.test(alpha_micro$Shannon, sample_data(micro_rare2)$Diet_type) #ns  chi-squared = 2.8402, df = 2, p-value = 0.2417
kruskal.test(div$diet_shannon, div$Species) #ns chi-squared=14.036, df=13, p=0.3713
kruskal.test(div$micro_shannon[which(!(div$Diet_type=='Intermediate'))], div$Diet_type[which(!(div$Diet_type=='Intermediate'))]) #still ns to compare generalist to specialist, chi-squared = 1.934, df = 1, p-value = 0.1643
kruskal.test(alpha_micro$Observed, sample_data(micro_rare2)$Diet_type) #ns

hist(div$micro_shannon[div$Diet_type=='Generalist'], breaks=10)
hist(div$micro_shannon[div$Diet_type=='Specialist'], breaks=10)
extremes<-div[which(!(div$Diet_type=='Intermediate')),] #54 individuals, just generalists and specailists
kruskal.test(extremes$micro_shannon, extremes$Diet_type) #ns chi-squared=1.934, df=1, p=0.1643, not significant even if intermediates are excluded

#pdf(file='figs_diet_and_microbiome2/diet_micro_scatter2017-2019.pdf', height=5, width=5)
par(mar=c(5, 4, 2, 6.5), xpd=TRUE)
plot(div$micro_shannon, div$diet_shannon, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (Shannon)', ylab='Diet alpha diversity (Shannon)',
     main='2017-2019')
abline(lm(div$diet_shannon~div$micro_shannon), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(div$micro_shannon, div$diet_shannon, method='kendall') #tau=0.02945736, p=0.6191

#mean alpha by species
par(mar=c(5, 4, 2, 6.5), xpd=TRUE)
plot(index_ab$alpha_micro, index_ab$alpha_diet, col=factor(rownames(index_ab)), lwd=2,
     xlab='Mean microbiome alpha diversity (Shannon)', ylab='Mean diet alpha diversity (Shannon)')
abline(lm(index_ab$alpha_diet~index_ab$alpha_micro), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(rownames(index_ab))), col=1:nlevels(factor(rownames(index_ab))), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(index_ab$alpha_micro, index_ab$alpha_diet, method='kendall')

pdf(file='figs_diet_and_microbiome2/unifrac_dietType2017-2019.pdf', height=5, width=6)
ordinate(micro_rare2, "PCoA", "unifrac") %>% 
  plot_ordination(micro_rare2, ., color='Diet_type', title = "UniFrac", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Diet type"))+
  scale_color_manual(values=c('#CC79A7', '#E69F00', '#56B4E9')) + geom_point(size=3) +
  theme_bw() + stat_ellipse(level=0.5, lty=2, lwd=0.6) +
  theme(axis.text = element_text(size = 11)) 
dev.off()

o<-ordinate(micro_rare2, "PCoA", "bray")
plot(o$vectors[1:130,1], -o$vectors[1:130,2])

pdf(file='figs_diet_and_microbiome2/bray_dietType2017-2019.pdf', height=5, width=5.5)
ordinate(micro_rare2, "PCoA", "bray") %>% 
  plot_ordination(micro_rare2, ., color='Diet_type', title = "", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Diet type"))+
  scale_color_manual(values=c('gray1', 'gray50', 'gray80'), labels=c("High diversity","Intermediate","Low diversity")) + geom_point(size=3) +
  theme_bw() + stat_ellipse(level=0.5, lty=2, lwd=0.6) +
  theme(axis.text = element_text(size = 11)) +
  scale_y_reverse()
dev.off()
#'#CC79A7', '#E69F00', '#56B4E9'
micro_rare2@sam_data$Year<-as.factor(micro_rare2@sam_data$Year) #make non-continuous

save.image('diet_and_microbiome2.RData')         

#mantel tests for phylosymbiosis comparison
identical(attr(bray_micro2, 'Labels'), attr(bray_diet2, 'Labels')) #yes!
mantel(bray_micro2, bray_diet2, method = "spearman", permutations = 999, na.rm = TRUE)
#Mantel statistic r: 0.08962 Significance: 0.044 
aa = as.vector(bray_micro2)
tt = as.vector(bray_diet2)
mat = data.frame(aa,tt)

pdf(file='figs_diet_and_microbiome2/mantel_dietBC2017-2019.pdf', height=5, width=5)
plot(mat$aa, mat$tt, col=alpha('black', 0.4), pch=16,  ylim=c(0,1), xlim=c(0,1),
     xlab='Microbiome similarity (Bray-Curtis)', ylab='Diet similarity (Bray-Curtis)', cex=1.3)
abline(lm(mat$tt~mat$aa), xpd=F, lty=2, lwd=1.5, col='red')
dev.off()

identical(attr(unifrac_micro, 'Labels'), attr(unifrac_diet, 'Labels')) #yes!
mantel(unifrac_micro, unifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.2248 , p=0.001
mat$bb = as.vector(unifrac_micro)
mat$cc = as.vector(unifrac_diet)

pdf(file='figs_diet_and_microbiome2/mantel_dietUni2017-2019.pdf', height=5, width=5)
plot(mat$bb, mat$cc, col=alpha('black', 0.4), pch=16, ylim=c(0.3,1), xlim=c(0.3,1),
     xlab='Microbiome similarity (UniFrac)', ylab='Diet similarity (UniFrac)', cex=1.3)
abline(lm(mat$cc~mat$bb), xpd=F, lty=2, lwd=1.5, col='red')
dev.off()

identical(attr(weightedu_micro, 'Labels'), attr(weightedu_diet, 'Labels')) #yes!
mantel(weightedu_micro, weightedu_diet, method = "spearman", permutations = 999, na.rm = TRUE) #-0.04771, p=0.777
mat$dd = as.vector(weightedu_micro)
mat$ee = as.vector(weightedu_diet)
plot(mat$dd, mat$ee, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (weighted UniFrac)', ylab='Diet similarity (weighted UniFrac)', cex=1.3)
abline(lm(mat$ee~mat$dd), xpd=F, lty=2, lwd=1.5, col='red')

identical(attr(jac_micro, 'Labels'), attr(jac_diet, 'Labels')) #yes!
mantel(jac_micro, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #0.21, p=0.001
mat$ff = as.vector(jac_micro)
mat$gg = as.vector(jac_diet)
plot(mat$ff, mat$gg, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (Jaccard)', ylab='Diet similarity (Jaccard)', cex=1.3)
abline(lm(mat$gg~mat$ff), xpd=F, lty=2, lwd=1.5, col='red')

#mantel tests for NY-subset
#distances for NY subset
bray_micro_ny<-distance(subset_samples(micro_rare2, State=='NY'), method="bray")
bray_diet_ny<-distance(subset_samples(diet_rare2, State=='NY'), method="bray")
jac_micro_ny<-distance(subset_samples(micro_rare2, State=='NY'), method="jaccard", binary=T)
jac_diet_ny<-distance(subset_samples(diet_rare2, State=='NY'), method="jaccard", binary=T)
uni_micro_ny<-distance(subset_samples(micro_rare2, State=='NY'), method="unifrac")
uni_diet_ny<-distance(subset_samples(diet_rare2, State=='NY'), method="unifrac")
wuni_micro_ny<-distance(subset_samples(micro_rare2, State=='NY'), method="wunifrac")
wuni_diet_ny<-distance(subset_samples(diet_rare2, State=='NY'), method="wunifrac")

identical(attr(bray_micro_ny, 'Labels'), attr(bray_diet_ny, 'Labels')) #yes!
mantel(bray_micro_ny, bray_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.1071, p=0.047
hh = as.vector(bray_micro_ny)
ii = as.vector(bray_diet_ny)
mat2 = data.frame(hh,ii)
plot(mat2$hh, mat2$ii, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (Bray-Curtis)', ylab='Diet similarity (Bray-Curtis)', cex=1.3)
abline(lm(mat2$ii~mat2$hh), xpd=F, lty=2, lwd=1.5, col='red')

identical(attr(jac_micro_ny, 'Labels'), attr(jac_diet_ny, 'Labels')) #yes!
mantel(jac_micro_ny, jac_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.2334, p=0.001
mat2$jj = as.vector(jac_micro_ny)
mat2$kk = as.vector(jac_diet_ny)
plot(mat2$jj, mat2$kk, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (Jaccard)', ylab='Diet similarity (Jaccard)', cex=1.3)
abline(lm(mat2$kk~mat2$jj), xpd=F, lty=2, lwd=1.5, col='red')

identical(attr(uni_micro_ny, 'Labels'), attr(uni_diet_ny, 'Labels')) #yes!
mantel(uni_micro_ny, uni_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) #0.2511, p=0.001
mat2$ll = as.vector(uni_micro_ny)
mat2$mm = as.vector(uni_diet_ny)
plot(mat2$ll, mat2$mm, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (UniFrac)', ylab='Diet similarity (UniFrac)', cex=1.3)
abline(lm(mat2$mm~mat2$ll), xpd=F, lty=2, lwd=1.5, col='red')

identical(attr(wuni_micro_ny, 'Labels'), attr(wuni_diet_ny, 'Labels')) #yes!
mantel(wuni_micro_ny, wuni_diet_ny, method = "spearman", permutations = 999, na.rm = TRUE) #-0.08638, p=0.893
mat2$nn = as.vector(wuni_micro_ny)
mat2$oo = as.vector(wuni_diet_ny)
plot(mat2$nn, mat2$oo, col=alpha('black', 0.4), pch=16,
     xlab='Microbiome similarity (weighted UniFrac)', ylab='Diet similarity (weighted UniFrac)', cex=1.3)
abline(lm(mat2$oo~mat2$nn), xpd=F, lty=2, lwd=1.5, col='red')

ordinate(micro_rare2, "PCoA", "bray") %>% 
  plot_ordination(micro_rare2, ., color='Diet_type', title = "", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Diet type"))+
  scale_color_manual(values=c('gray1', 'gray50', 'gray80'), labels=c("High diversity","Intermediate","Low diversity")) + geom_point(size=3) +
  theme_bw() + stat_ellipse(level=0.95, lty=2, lwd=0.6) +
  theme(axis.text = element_text(size = 11)) +
  scale_y_reverse()


##diet exluding 4 species that differ between batch 1 and full
include<-c('MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','CAWA','NAWA','OVEN')
kruskal.test(alpha_micro$Shannon[which(sample_data(micro_rare2)$Species %in% include)], 
             sample_data(micro_rare2)$Diet_type[sample_data(micro_rare2)$Species %in% include]) #ns

#diet index with other new alpha and beta metrics
#original calculate indexes
#i1 is shannon+BC
#i2 is shannon+uni
#calculate chao1 and faith pd for diet alpha
alpha_diet$chao<-estimate_richness(diet_rare2, measures=c('Chao1'))[1] #add to alpha dataframe
library(btools) #for faiths pd
testpd<-estimate_pd(diet_rare2)
alpha_diet$PD<-testpd$PD
#correlation between alphas
plot(alpha_diet$Shannon~alpha_diet$chao$Chao1)
plot(alpha_diet$Shannon~alpha_diet$PD)
plot(alpha_diet$chao$Chao1~alpha_diet$PD)
#for supplemental table
cor.test(alpha_diet$Shannon, alpha_diet$PD) #pearson r=0.704, p<0.001
cor.test(alpha_diet$Shannon, alpha_diet$chao$Chao1) #r=0.394, p<0.001
cor.test(alpha_diet$PD, alpha_diet$chao$Chao1) #r=0.708, p<0.001

pdf(file='figs_diet_and_microbiome2/diet_index_correlations/shan_pd.pdf', height=5, width=5)
plot(alpha_diet$Shannon, alpha_diet$PD, pch=16,
     xlab='Shannon index (diet)' , ylab='Faith\'s pd (diet)', cex=1.3)
abline(lm(alpha_diet$PD~alpha_diet$Shannon), lty=2, lwd=1.5, col='gray')
dev.off()
pdf(file='figs_diet_and_microbiome2/diet_index_correlations/shan_chao.pdf', height=5, width=5)
plot(alpha_diet$Shannon, alpha_diet$chao$Chao1, pch=16,
     xlab='Shannon index (diet)' , ylab='Chao1 (diet)', cex=1.3)
abline(lm(alpha_diet$chao$Chao1~alpha_diet$Shannon), lty=2, lwd=1.5, col='gray')
dev.off()
pdf(file='figs_diet_and_microbiome2/diet_index_correlations/pd_chao.pdf', height=5, width=5)
plot(alpha_diet$PD, alpha_diet$chao$Chao1, pch=16,
     xlab='Faith\'s PD (diet)' , ylab='Chao1 (diet)', cex=1.3)
abline(lm(alpha_diet$chao$Chao1~alpha_diet$PD), lty=2, lwd=1.5, col='gray')
dev.off()

#Jaccard distance
jac_diet<-distance(diet_rare2, method='jaccard', binary=T)
jac_micro<-distance(micro_rare2, method='jaccard', binary=T)


mantel(bray_diet2, unifrac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.4912, p=0.001
mantel(bray_diet2, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.6574, p=0.001
mantel(unifrac_diet, jac_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.7306, p=0.001
mantel(bray_diet2, weightedu_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.3771, p=0.001
mantel(jac_diet, weightedu_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.1992, p=0.001
mantel(unifrac_diet, weightedu_diet, method = "spearman", permutations = 999, na.rm = TRUE) #r=0.1289, p=0.001

##correlation between alpha diet and alpha micro w/pd and chao
div$diet_PD<-testpd$PD
div$diet_chao1<-alpha_diet$chao$Chao1
#calculate micro pd and chao
alpha_micro$chao<-estimate_richness(micro_rare2, measures=c('Chao1'))[1] #add to alpha dataframe
testpd2<-estimate_pd(micro_rare2)
alpha_micro$PD<-testpd2$PD
div$micro_PD<-testpd2$PD
div$micro_chao1<-alpha_micro$chao$Chao1
plot(div$micro_PD, div$diet_PD, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (PD)', ylab='Diet alpha diversity (PD)',
     main='2017-2019')
abline(lm(div$diet_PD~div$micro_PD), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(div$micro_PD, div$diet_PD, method='kendall') #tau=0.1096005, p=0.06438
plot(div$micro_chao1, div$diet_chao1, col=factor(div$Species), lwd=2,
     xlab='Microbiome alpha diversity (Chao1)', ylab='Diet alpha diversity (Chao1)',
     main='2017-2019')
abline(lm(div$diet_chao1~div$micro_chao1), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=1, pt.lwd=2, inset=c(-0.35,0))
dev.off()
cor.test(div$micro_chao1, div$diet_chao1, method='kendall') #tau=0.1232518, p=0.03797

#adonis on diet type
adonis2(jac_micro~Diet_type, data=md) #r2=0.0197, p=0.001***
permutest(betadisper(jac_micro, md$Diet_type)) # F=9.0666, p=0.001**
plot(betadisper(jac_micro, md$Diet_type))
pal<-palette(c('gray1', 'gray50', 'gray80'))

myplotbetadisper <- function (x, axes = c(1, 2), cex = 0.7, pch = seq_len(ng), col = NULL, 
                              lty = "solid", lwd = 1, hull = TRUE, ellipse = FALSE, conf, 
                              segments = TRUE, seg.col = "grey", seg.lty = lty, seg.lwd = lwd, 
                              label = TRUE, label.cex = 1, ylab, xlab, main, sub, 
                              fillrect="white", coltextrect="black", alphaPoints=0.5, 
                              labPoints=NULL, poslabPoints=4, ...) 
{
  localAxis <- function(..., col, bg, pch, cex, lty, lwd) axis(...)
  localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
  localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
  Ellipse <- function(scrs, centres, conf, col, lty, lwd, ...) {
    mat <- cov.wt(scrs, center = centres)
    if (mat$n.obs == 1) 
      mat$cov[] <- 0
    xy <- if (mat$n.obs > 1) {
      vegan:::veganCovEllipse(mat$cov, mat$center, conf)
    }
    else {
      scrs
    }
    vegan:::ordiArgAbsorber(xy, FUN = lines, col = col, lty = lty, 
                            lwd = lwd, ...)
  }
  if (missing(main)) 
    main <- deparse(substitute(x))
  if (missing(sub)) 
    sub <- paste("method = \"", attr(x, "method"), "\"", 
                 sep = "")
  if (missing(xlab)) 
    xlab <- paste("PCoA", axes[1])
  if (missing(ylab)) 
    ylab <- paste("PCoA", axes[2])
  t <- if (missing(conf)) {
    1
  }
  else {
    sqrt(qchisq(conf, df = 2))
  }
  g <- scores(x, choices = axes)
  ng <- length(levels(x$group))
  lev <- levels(x$group)
  if (is.null(col)) {
    col <- palette()
  }
  col <- rep_len(col, ng)
  colpts <- apply(col2rgb(col), 2, addAlpha, alpha=alphaPoints)
  seg.col <- rep_len(seg.col, ng)
  plot(g$sites, asp = 1, type = "n", axes = FALSE, ann = FALSE, 
       ...)
  if (is.matrix(g$centroids)) {
    for (i in seq_along(lev)) {
      curlev <- lev[i]
      take <- x$group == curlev
      j <- which(lev == curlev)
      if (segments) {
        segments(g$centroids[j, 1L], g$centroids[j, 2L], 
                 g$sites[take, 1L], g$sites[take, 2L], col = seg.col[i], 
                 lty = seg.lty, lwd = seg.lwd)
      }
      if (hull) {
        ch <- chull(g$sites[take, ])
        ch <- c(ch, ch[1])
        lines(x$vectors[take, axes][ch, ], col = col[i], 
              lty = lty, lwd = lwd, ...)
      }
      if (ellipse) {
        Ellipse(g$sites[take, , drop = FALSE], centres = g$centroids[j, 
        ], conf = t, col = col[i], lty = lty, lwd = lwd, 
        ...)
      }
      points(g$centroids[j, , drop = FALSE], pch = 16, 
             cex = 1, col = col[i], ...)
    }
  }
  else {
    if (segments) {
      segments(g$centroids[, 1L], g$centroids[, 2L], g$sites[, 
                                                             1L], g$sites[, 2L], col = seg.col, lty = seg.lty, 
               ...)
    }
    if (hull) {
      ch <- chull(g$sites)
      ch <- c(ch, ch[1])
      lines(x$vectors[, axes][ch, ], col = col[1L], lty = lty, 
            lwd = lwd, ...)
    }
    if (ellipse) {
      Ellipse(g$sites, centres = g$centroids, conf = t, 
              col = col[1L], lty = lty, lwd = lwd, ...)
    }
    points(g$centroids[, 1L], g$centroids[, 2L], pch = 16, 
           cex = 1, col = col[1L], ...)
  }
  points(g$sites, pch = pch[x$group], cex = cex, col = col[x$group], 
         ...)
  if (!is.null(labPoints)) {
    text(g$sites, labels=labPoints, pos=poslabPoints,
         cex = cex, col = col[x$group])
  }
  if (label) {
    myordilabel(x, labels=c('High diversity','Intermediate diversity','Low diversity' ), display = "centroids", choices = axes, cex = label.cex, fill=fillrect, col=coltextrect)
  }
  localTitle(main = main, xlab = xlab, ylab = ylab, sub = sub, 
             ...)
  localAxis(1, ...)
  localAxis(2, ...)
  localBox(...)
  class(g) <- "ordiplot"
  invisible(g)
}

pdf(file='figs_diet_and_microbiome2/betadisper_jac_diettype.pdf', height=5, width=5)
myplotbetadisper(betadisper(jac_micro, md$Diet_type), ellipse = TRUE, hull = FALSE, 
                 fillrect=col.fill.rect, seg.lwd = 0.5, cex=0.8, lwd=1.4,
                 alphaPoints=0.7, labPoints=NULL, sub='',
                 main= "Jaccard", coltextrect=c('black','gray','lightgray'), label.cex=0.9)
dev.off()
pdf(file='figs_diet_and_microbiome2/betadisper_unifrac_diettype.pdf', height=5, width=5)
myplotbetadisper(betadisper(unifrac_micro, md$Diet_type), ellipse = TRUE, hull = FALSE, 
                 fillrect=col.fill.rect, seg.lwd = 0.5, cex=0.8, lwd=1.4,
                 alphaPoints=0.7, labPoints=NULL, sub='',
                 main= "UniFrac", coltextrect=c('black','gray','lightgray'), label.cex=0.9)
dev.off()

#plot of diet similarity 
pal<-palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
               "red2", "#75DED6", "#75E563", "black", "#B44AE1"))
pdf(file="~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/figs/pcoa_diet.pdf", height=5, width=5)
ordinate(diet_rare2, "PCoA", "bray") %>% 
  plot_ordination(micro_rare2, ., color='Species', shape = "State", title = "", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Species"))+
  scale_color_manual(values=pal) + geom_point(size=3) +
  theme_bw() +  labs(shape="Locality", colour="Species") +
  theme(axis.text = element_text(size = 11))
dev.off()

#save individual diet matricies to file
saveRDS(bray_diet2, "bray_diet2.rds")
saveRDS(jac_diet, 'jac_diet.rds')
saveRDS(unifrac_diet, "unifrac_diet.rds")
saveRDS(weightedu_diet, "weightedu_diet.rds")
saveRDS(jac_diet, "jac_diet.rds")
saveRDS(bray_micro2, 'bray_micro2.rds')
saveRDS(jac_micro, 'jac_micro.rds')
saveRDS(unifrac_micro, "unifrac_micro.rds")
saveRDS(weightedu_micro, 'weightedu_micro.rds')

#kruskal tests of microbiome by diet type 
kruskal.test(alpha_micro$chao$Chao1, sample_data(micro_rare2)$Diet_type) #ns  chi-squared = 5.3879, df = 2, p-value = 0.06761
kruskal.test(alpha_micro$PD, sample_data(micro_rare2)$Diet_type) #ns  chi-squared = 4.4195, df = 2, p-value = 0.1097

kruskal.test(div$micro_chao1[which(!(div$Diet_type=='Intermediate'))], div$Diet_type[which(!(div$Diet_type=='Intermediate'))]) #still ns to compare generalist to specialist, chi-squared = 1.7914, df = 1, p-value = 0.1808
kruskal.test(div$micro_PD[which(!(div$Diet_type=='Intermediate'))], div$Diet_type[which(!(div$Diet_type=='Intermediate'))]) #still ns to compare generalist to specialist, chi-squared = 1.2731, df = 1, p-value = 0.2592

#make scatter of chao1 for figure 2b
pdf(file='figs_diet_and_microbiome2/diet_micro_scatter2017-2019_chao1.pdf', height=5, width=5)
par(mar=c(5, 4, 2, 7.3), xpd=TRUE)
pal<-palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
               "red2", "#75DED6", "#75E563", "black", "#B44AE1"))
plot(div$micro_chao1, div$diet_chao1, col=factor(div$Species), pch=16,
     xlab='Microbiome alpha diversity (Chao1)', ylab='Diet alpha diversity (Chao1)', cex=1.3)
abline(lm(div$diet_chao1~div$micro_chao1), xpd=F, lty=2, lwd=1.5)
legend("right", legend=levels(factor(div$Species)), col=1:nlevels(factor(div$Species)), pch=16, inset=c(-0.37,0), pt.cex=1.2,bty='n')
dev.off()


saveRDS(bray_micro_ny, 'bray_micro_ny.rds')
saveRDS(bray_diet_ny, 'bray_diet_ny.rds')
saveRDS(jac_micro_ny, 'jac_micro_ny.rds')
saveRDS(jac_diet_ny, 'jac_diet_ny.rds')
saveRDS(uni_micro_ny, 'uni_micro_ny.rds')
saveRDS(uni_diet_ny, 'uni_diet_ny.rds')
saveRDS(wuni_micro_ny, 'wuni_micro_ny')
saveRDS(wuni_diet_ny, 'wuni_diet_ny.rds')

#PGLS
wtre<-read.tree("~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/eliot_example/scaledWarbler_full.tre")
wtre$tip.label<-c('OVEN','WEWA','BAWW','NAWA','COYE','HOWA','AMRE','BTBW','NOPA','MAWA','BLBW','CSWA','MYWA','BTNW','CAWA')
wtre$root.edge<-0
wtre2<-drop.tip(wtre, 'WEWA')
wtre2$node.label<-NULL
library(phylolm)
#cor.test was (div$micro_chao1, div$diet_chao1)

fit_shan = phylolm(alpha_micro~alpha_diet,data=index_ab,phy=wtre,model="lambda", boot=1000) 
summary(fit_shan)

#add species averages for other alpha metrics
index_ab$diet_chao1=c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14)
index_ab$micro_chao1=c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)
index_ab$diet_pd=c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14)
index_ab$micro_pd=c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14)

fit_chao = phylolm(micro_chao1~diet_chao1,data=index_ab,phy=wtre,model="lambda", boot=1000) 
summary(fit_chao)

fit_pd = phylolm(micro_pd~diet_pd,data=index_ab,phy=wtre,model="lambda", boot=1000) 
summary(fit_pd)

#standard deviation for N individuals per species in these analyses
micro_rare2
sd(table(micro_rare2@sam_data$Species))

#make new diet type plots with species in color
#add diet type to div dataframe
identical(rownames(div), rownames(micro_rare2@sam_data))
div$Diet_type<-micro_rare2@sam_data$Diet_type
pdf(file='figs_diet_and_microbiome2/boxplot_shannon_diettype2017-2019_color.pdf', height=3.5, width=3.5)
div %>% ggplot(aes(x=Diet_type, y=micro_shannon))+  
  geom_boxplot(show.legend = F, aes(fill=Diet_type), outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0), 
             pch=21, aes(fill=Species), show.legend = F, size=2.5) +
  scale_x_discrete(limits = rev(levels(factor(micro_rare2@sam_data$Diet_type))),labels=c("Low\ndiversity","Intermediate","High\ndiversity")) +
  scale_fill_manual(values = c("Specialist" = "gray80","Intermediate" = "gray50","Generalist" = "gray30",
                               'AMRE'='#6D7BD0','BAWW'='blue','BLBW'='#D8C278','BTBW'='#D5E356','BTNW'='#D45691',
                               'CAWA'='#D988DF','COYE'='#D4E1CE','CSWA'='#79ACCD','HOWA'='#C1958E',
                               'MAWA'='red2','MYWA'='#75DED6','NAWA'='#75E563','NOPA'='black','OVEN'='#B44AE1')) +
  xlab('Diet type') + ylab('Microbiome alpha diversity (Shannon)') +
  ylim(0,6.5) + theme_bw() +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.text.x = element_text(size = 11))
dev.off()

#and make new PCOA diet type plot w/species by color
pdf(file='figs_diet_and_microbiome2/bray_dietType2017-2019_color.pdf', height=5, width=5.5)
ordinate(micro_rare2, "PCoA", "bray") %>% 
  plot_ordination(micro_rare2, ., color='Species', shape = "Diet_type", title = "", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Species"))+
  scale_color_manual(values = c("Specialist" = "gray80","Intermediate" = "gray50","Generalist" = "gray1",
                                'AMRE'='#6D7BD0','BAWW'='blue','BLBW'='#D8C278','BTBW'='#D5E356','BTNW'='#D45691',
                                'CAWA'='#D988DF','COYE'='#D4E1CE','CSWA'='#79ACCD','HOWA'='#C1958E',
                                'MAWA'='red2','MYWA'='#75DED6','NAWA'='#75E563','NOPA'='black','OVEN'='#B44AE1')) + 
  geom_point(size=3) +
  theme_bw() +  labs(shape="Diet type", colour="Species") +
  theme(axis.text = element_text(size = 11)) + 
  stat_ellipse(geom = 'path' , aes(color = Diet_type, group=Diet_type,linetype=factor(micro_rare2@sam_data$Diet_type)), level=0.5, lwd=0.6) +
  scale_y_reverse() + scale_linetype_manual(values=c(2,2,2), guide='none') 
dev.off()  
