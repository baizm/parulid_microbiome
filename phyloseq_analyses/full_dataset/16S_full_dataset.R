setwd("~/Documents/Toews_Lab/16s_merged/")
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
library(lubridate)
library(VennDiagram)

#read in merged data from qiime2
########-------create physeq object---------------------########
####------read in OTU table, has to be a matrix---------------------
otu_table<-read.csv('qiime_output/feature-table.tsv', sep='\t', header=T, skip=1, check.names=FALSE)
#rearrange columns in abc order (don't know if necessary)
otu_table<-otu_table[,order(colnames(otu_table))]
#convert to matrix
otumat<-as.matrix(otu_table[,2:408])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
###-------read in taxonomy table-----------------
tax_table<-read.csv('qiime_output/silva_taxonomy_merged.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('metadata_16s_merged.csv', header=T,)
#put rows in abc order
sampledata<-sampledata[order(sampledata$id),]
rownames(sampledata)<-sampledata$id
x<-sampledata$id[(which(!(sampledata$id %in% colnames(otumat))))]
sampledata<-sampledata[which(!(sampledata$id %in% x)),] #get rid of samples filtered by qiime (those in x)
sampledata<-sampledata[,2:31] #get rid of id column
###-------read in nwk tree with ape package-----------
tree<-read.tree('qiime_output/tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
physeq = phyloseq(OTU, TAX, SAM, tree)
#sums of OTUs across samples
taxa_sums(physeq)
plot_richness(physeq, x="Species", color="Species",measures="Shannon") + theme(legend.position="none") 

table(sampledata$Species)

#get rid of negatives, oven experiment, and vermivora, and recaps
negs<-rownames(sampledata[which(!(sampledata$Type=='sample')),]) ###negatives to exclude
ovex<-rownames(sampledata[which(sampledata$Method=='Frozen'),]) ###experimental ovens to exclude
v<-c('BRWA','BWWA','GWWA')
verm<-rownames(sampledata[which(sampledata$Species %in% v),]) ###vermivora to exclude
recaps<-sampledata[which(sampledata$Band %in% unique(sampledata$Band[duplicated(sampledata$Band)])),] #identify recaps
recaps<-recaps[6:35,] #remove negatives
recaps<-recaps[which(!(recaps$Species %in% v)),] #remove vermivora
recaps<-recaps[which(!(recaps$Method=='Frozen')),] #remove experimental ovens
recaps<-recaps[which(!(recaps$Species=='OVEN')),] #remove other OVENs, no longer duplicates
rows <- sample(nrow(recaps)) #shuffle row indices
recaps2<-recaps[rows,] #use random vector to reorder recaps
recaps3<-recaps2[which(duplicated(recaps2$Band)),]
recaps4<-rownames(recaps3) ###recaps to remove

remove<-c(negs,ovex,verm,recaps4)
physeq2<-prune_samples(rownames(sampledata[which(!(rownames(sampledata) %in% remove)),]), physeq)

p17<-prune_samples(rownames(sampledata[which(sampledata$Year=='2017'),]), physeq2)
p17<-prune_taxa(taxa_sums(p17) > 0, p17)
p18<-prune_samples(rownames(sampledata[which(sampledata$Year=='2018'),]), physeq2)
p18<-prune_taxa(taxa_sums(p18) > 0, p18)
p19<-prune_samples(rownames(sampledata[which(sampledata$Year=='2019'),]), physeq2)
p19<-prune_taxa(taxa_sums(p19) > 0, p19)
p20<-prune_samples(rownames(sampledata[which(sampledata$Year=='2020'),]), physeq2)
p20<-prune_taxa(taxa_sums(p20) > 0, p20)

table(sampledata$Year)
#overlap of OTUs by year
venn.diagram(
  x = list(rownames(p17@otu_table),rownames(p20@otu_table),rownames(p18@otu_table),rownames(p19@otu_table)),
  category.names = c("2017" , "2020", '2018', '2019'),
  filename = 'otu_venn_by_year.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c('indianred1','darkolivegreen','tan4','darkblue')
) #221 OTUs overlap in all 4 years

#overlapping otu abundance in 2017
o17<-p17@otu_table[which(rownames(p17@otu_table) %in% overlap_otus),]
os17<-taxa_sums(o17)
sum(os17)/sum(p17@otu_table) #84% rel abundance are overlapping OTUs

o18<-p18@otu_table[which(rownames(p18@otu_table) %in% overlap_otus),]
os18<-taxa_sums(o18)
sum(os18)/sum(p18@otu_table) #78% rel abundance are overlapping OTUs

o19<-p19@otu_table[which(rownames(p19@otu_table) %in% overlap_otus),]
os19<-taxa_sums(o19)
sum(os19)/sum(p19@otu_table) #86% rel abundance are overlapping OTUs

o20<-p20@otu_table[which(rownames(p20@otu_table) %in% overlap_otus),]
os20<-taxa_sums(o20)
sum(os20)/sum(p20@otu_table) #68% rel abundance are overlapping OTUs

df <- data.frame(year=rep(c('2017', '2018', '2019', '2020'), each=2),
                 otus=rep(c('Overlapping (1164 otus)', 'Non-overlapping'), times=4),
                 rel.abundance=c(sum(os17)/sum(p17@otu_table), 1-sum(os17)/sum(p17@otu_table),
                         sum(os18)/sum(p18@otu_table), 1-sum(os18)/sum(p18@otu_table),
                         sum(os19)/sum(p19@otu_table), 1-sum(os19)/sum(p19@otu_table),
                         sum(os20)/sum(p20@otu_table), 1-sum(os20)/sum(p20@otu_table)))
pdf('rel.abundance_overlappingOTUs.pdf', height=5, width=4)
ggplot(df, aes(fill=otus, y=rel.abundance, x=year)) + 
  geom_bar(position='stack', stat='identity') +
  scale_fill_manual('rel.abundance', values=c('darkgray', 'steelblue'))
dev.off() 


#do rarecurve
#decide which species to remove from larger analyses (bc <5 individuals)

####---------find samples to remove due to low read counts-------------
#make data frame with total read counts after filtering using sample_sums
sdt = data.table(as(sample_data(physeq2), "data.frame"),
                 TotalReads = sample_sums(physeq2), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
summary(sdt$TotalReads)
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
pdf(file='figs_for_paper/readDepth_byLane.pdf', height=5, width=4)
boxplot(sdt$TotalReads~sdt$Sequencing_lane, outline=F, xlab='Sequencing run', ylab='N reads/sample', boxwex=0.5, col='gray94')
stripchart(sdt$TotalReads[which(sdt$Sequencing_lane=='1')],vertical=T,pch=21, method='jitter',add=T,bg=alpha('gray55', alpha=0.6))
stripchart(sdt$TotalReads[which(sdt$Sequencing_lane=='2')],vertical=T,pch=21, method='jitter',add=T,bg=alpha('gray18', alpha=0.6),at=2)
dev.off()
sum(sdt$TotalReads[which(sdt$Sequencing_lane=='1')]) #3,817,515
sum(sdt$TotalReads[which(sdt$Sequencing_lane=='1')])/length(sdt$TotalReads[which(sdt$Sequencing_lane=='1')]) #15645.55
sum(sdt$TotalReads[which(sdt$Sequencing_lane=='2')]) #2,644,775
sum(sdt$TotalReads[which(sdt$Sequencing_lane=='2')])/length(sdt$TotalReads[which(sdt$Sequencing_lane=='2')]) #22413.35
#look at new rarecurve
#pdf(file='rarecurve_merged.pdf', height=4, width=4)
rarecurve(t(otu_table(physeq2)), step=50, lwd=0.8, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,9000))
dev.off()
#data frame of samples with reads below threshold
high<-sdt[which(sdt$TotalReads>4000),]
barplot(table(physeq2@sam_data$Species), las=2, ylim=c(0,40))
barplot(table(high$Species), las=2, ylim=c(0,40))
low<-sdt[which(sdt$TotalReads<4000),]
table(low$Species) 
table(physeq2@sam_data$Species) #for comparison, it looks like mostly removes from species that already have a lot of individuals
dim(low) 

#make vector of samples to keep
keep_samples<-sdt$SampleID[which(sdt$SampleID %in% high$SampleID)]
#make new phyloseq object to exclude low read samples
physeq3<-prune_samples(keep_samples, physeq2)
#plot original vs pruned samples by species
#pdf(file='samples_by_species.pdf', height=8, width=5)
par(mfrow=c(2,1), mar=c(4,4,2,1))
barplot(table(physeq2@sam_data$Species), las=2, ylab='N', main='Original', ylim=c(0,40))
abline(h=5, col='red', lty=2)
barplot(table(physeq3@sam_data$Species), las=2, ylab='N', main='Pruned', ylim=c(0,40))
abline(h=5, col='red', lty=2)
#dev.off()
min(sample_sums(physeq3)) 
max(sample_sums(physeq3)) 

#write file of pruned warblers for qiime filtering
#prune species with < 5 indiv 
l4<-c('BLPW','LOWA','NOWA','PIWA','YEWA','YPWA') #species with less than 5 indiv
#new pruned dataset removing poorly sampled species
physeq4<-prune_samples(!(grepl(paste(l4,collapse='|'),sample_data(physeq3)$Species)), physeq3)
#write.csv(rownames(physeq4@sam_data), './qiime_diversity/samples_to_keep.txt',quote=F,row.names=F)

###--------identify the core microbiome, general warbler GM -----------
core<-data.frame(physeq4@otu_table) #make a copy of non-rarefied OTU table
colnames(core)<-colnames(physeq4@otu_table)
core[core==0]<-NA #turn zeros to NA
core$sumNA<-rowSums(is.na(core)) #add column of sums of NAs per row
summary(core$sumNA) #min is 104, meaning present in (270-104=) 166/270, or 61% of individuals
core2<-core[which(core$sumNA<135),] #in more than half of the individuals, only 3 OTUs: 63afe8e6aac58bf0d670a82ca5bc574c,cc761daf51f27c423da57f3f1f0ff5cc,556864a5da3a811b67be9fc73488e926
core3<-core[order(core$sumNA),]
core3$prev<-c((270-core3$sumNA)/270) #add column of prevalence (proportion of individuals with this ASV)
core3$sum<-rowSums(core3[,1:270], na.rm=T)
which(core3$sumNA==270)
zo1<-which(core3$sum==0) #same as above. indexes of rows, remove all these OTUs
ro1<-rownames(core3[zo1,]) #otus to remove from physeq4 bc prevalence is zero

pdf(file='figs_for_paper/prevalence1.pdf', height=5, width=4)
hist(core3$prev,breaks=50, xlab="prevalence across individuals",main='')
dev.off()
pdf(file='figs_for_paper/prevalence2.pdf', height=4, width=4)
hist(core3$prev,breaks=50, xlim=c(0.1,0.65), ylim=c(0,30),
     xlab="prevalence across individuals",main='') #zoomed in
dev.off()

core3[1:10,270:273]
physeq4@tax_table[which(rownames(physeq4@tax_table)=='63afe8e6aac58bf0d670a82ca5bc574c'),]
physeq4@tax_table[which(rownames(physeq4@tax_table)=='556864a5da3a811b67be9fc73488e926'),]
physeq4@tax_table[which(rownames(physeq4@tax_table)=='cc761daf51f27c423da57f3f1f0ff5cc'),]
physeq4@tax_table[which(rownames(physeq4@tax_table)=='e09f13f193b7904630a9777d5428c247'),]

##----------rarefy------------------
sum(ps@otu_table)
pdf(file='sampling_byYear.pdf', height=5, width=4)
barplot(table(ps@sam_data$Year), ylab='N invidiuals', ylim=c(0,120))
dev.off()

ps<-rarefy_even_depth(physeq4,rngseed = 999, replace = F, trimOTUs = T, verbose = T) #270 individuals

#basic stats:
prune_taxa(taxa_sums(ps) > 0, ps) #note: checked that it does not include zero sum taxa!
dim(ps@tax_table) #12,048 ASVs
table(ps@tax_table[,2]) 
dim(table(ps@tax_table[,2])) #39 phyla
table(ps@sam_data$Species) #15 species
mean(table(ps@sam_data$Species)) #mean 18 individuals per species

# relative abundance top 8 phyla
topm<-as.data.frame(ps@tax_table)
sort(table(topm$Phylum), decreasing=T)
topmm<-names(sort(table(topm$Phylum), decreasing=T)[1:8]) #top 8 orders

ps2<-tax_glom(ps, taxrank = "Phylum", NArm=F) #glom by phylum
mps2<-psmelt(ps2)
mps2$Phylum[which(!(mps2$Phylum %in% topmm))]<-'other' #change non-common phyla to other 
table(mps2$Phylum)

#find spp which have top ASV
table(mps2$sample_Species[mps2$OTU=='63afe8e6aac58bf0d670a82ca5bc574c']) #represented in all 15 spp.

sum(mps2$Abundance) #1087560
sum(mps2$Abundance[which(mps2$Phylum=='p__Proteobacteria')]) #654871 out of 1087560 : 60%
sum(mps2$Abundance[which(mps2$Phylum=='p__Planctomycetota')]) #43890 out of 1087560 : 4%
sum(mps2$Abundance[which(mps2$Phylum=='p__Armatimonadota')]) #29088 out of 1087560 : 2.7%
sum(mps2$Abundance[which(mps2$Phylum=='p__Actinobacteriota')]) #70789 out of 1087560 : 6.5%
sum(mps2$Abundance[which(mps2$Phylum=='p__Bacteroidota')]) #25679 out of 1087560 : 2.3%
sum(mps2$Abundance[which(mps2$Phylum=='p__Verrucomicrobiota')]) #29090 out of 1087560 : 2.7%
sum(mps2$Abundance[which(mps2$Phylum=='p__Acidobacteriota')]) #29643 out of 1087560 : 2.7%
sum(mps2$Abundance[which(mps2$Phylum=='p__Firmicutes')]) #144814 out of 1087560 : 13%
#top 3 phyla 654871+70789+144814/1087560

#find abundance of Clostridia -- most abundant in Skeen paper
ps3<-tax_glom(ps, taxrank = "Class", NArm=F) #glom by class
mps3<-psmelt(ps3)
sum(mps3$Abundance)
sum(mps3$Abundance[which(mps3$Class=='c__Clostridia')]) #18608 out of 1087560 :1.7% 

mps2$Phylum2<-as.factor(mps2$Phylum)
mps2$Phylum2<-relevel(mps2$Phylum2, 'other')
pdf(file='figs_for_paper/abundance_microbiome_2017-2020.pdf', height=5, width=4)
level_order<-c('BTNW','MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','HOWA','CAWA','COYE','NAWA','BAWW','WEWA','OVEN')
ggplot(mps2, aes(x = factor(sample_Species, level=level_order), y = Abundance, fill = Phylum2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','turquoise4',
                             'rosybrown3','turquoise2','thistle1',
                             'rosybrown2','wheat4','sandybrown','black')) +
  ggtitle('') +
  scale_y_continuous(expand = c(0,0))  +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(axis.title.x=element_blank()) +
  scale_x_discrete(expand = c(0,0)) + theme(axis.text = element_text(size = 11)) +
  theme(legend.position="none")
dev.off()


#calculate alpha diversity
plot_richness(ps, x='Species', color='Species',measures=c('Observed','Shannon')) + 
  geom_boxplot() +
  theme(legend.position="none")
plot_richness(ps, x='State', color='State',measures=c('Observed','Shannon')) + 
  geom_boxplot() +
  theme(legend.position="none")
#pdf(file='alpha_byLane.pdf', height=5, width=4)
plot_richness(ps, x='Sequencing_lane', color='Sequencing_lane',measures=c('Observed','Shannon')) + 
  geom_boxplot() +
  theme(legend.position="none")
dev.off()

alpha<-estimate_richness(ps, split=T, measures=c('Observed', 'Shannon'))
rownames(alpha)<-rownames(ps@sam_data)
kruskal.test(alpha$Observed, sample_data(ps)$Species) 
kruskal.test(alpha$Shannon, sample_data(ps)$Species) 

kruskal.test(alpha$Observed, sample_data(ps)$Sequencing_lane)
kruskal.test(alpha$Shannon, sample_data(ps)$Sequencing_lane)
kruskal.test(alpha$Observed, sample_data(ps)$State)
kruskal.test(alpha$Observed, sample_data(ps)$Year)

#calculate distance matrices
#calculate distance matrix non-rare
bray<-distance(ps, method='bray')
#and unifrac
uni<-distance(ps, method='unifrac')
#and weighted unifrac
wuni<-distance(ps, method='wunifrac')

#pdf(file='beta_bySpecies.pdf', height=5, width=5)
ordinate(ps, "PCoA", "bray") %>% 
  plot_ordination(ps, ., color='Species', title = "Bray-Curtis (rarefied)", type='samples') +
  scale_color_manual(values=c('antiquewhite4', 'aquamarine4', 'azure4','black','blue4',
                              'brown1','burlywood','chartreuse4','cornflowerblue','darkgoldenrod',
                              'darkgray','darkviolet','khaki','lightpink4','darkred')) +
  geom_point(size=2) 
ordinate(ps, "PCoA", "uni") %>% 
  plot_ordination(ps, ., color='Species', title = "Unifrac (rarefied)", type='samples') +
  scale_color_manual(values=c('antiquewhite4', 'aquamarine4', 'azure4','black','blue4',
                              'brown1','burlywood','chartreuse4','cornflowerblue','darkgoldenrod',
                              'darkgray','darkviolet','khaki','lightpink4','darkred')) +
  geom_point(size=2) 
ordinate(ps, "PCoA", "wuni") %>% 
  plot_ordination(ps, ., color='Species', title = "Weighted Unifrac (rarefied)", type='samples') +
  scale_color_manual(values=c('antiquewhite4', 'aquamarine4', 'azure4','black','blue4',
                              'brown1','burlywood','chartreuse4','cornflowerblue','darkgoldenrod',
                              'darkgray','darkviolet','khaki','lightpink4','darkred')) +
  geom_point(size=2) 
dev.off()  

#for PERMANOVAS
md<-data.frame(ps@sam_data, row.names=rownames(ps@sam_data))
md$Year<-as.factor(md$Year) #make non-continuous
ps@sam_data$Year<-as.factor(ps@sam_data$Year) #make non-continuous
md$Sequencing_lane<-as.factor(md$Sequencing_lane) #make non-continuous
ps@sam_data$Sequencing_lane<-as.factor(ps@sam_data$Sequencing_lane) #make non-continuous

#pdf(file='beta_byLane.pdf', height=5, width=5)
par(mfrow=c(1,3), mar=c(4,4,2,1))
ordinate(ps, "PCoA", "bray") %>% 
  plot_ordination(ps, ., color='Sequencing_lane', title = "Bray-Curtis (rarefied)", type='samples') +
  geom_point(size=2) 
ordinate(ps, "PCoA", "uni") %>% 
  plot_ordination(ps, ., color='Sequencing_lane', title = "Unifrac (rarefied)", type='samples') +
  geom_point(size=2) 
ordinate(ps, "PCoA", "wuni") %>% 
  plot_ordination(ps, ., color='Sequencing_lane', title = "Weighted Unifrac (rarefied)", type='samples') +
  geom_point(size=2) 
dev.off()

pdf(file='figs_for_paper/beta_byLane.pdf', height=5, width=5)
par(mfrow=c(1,3), mar=c(4,4,2,1))
ordinate(ps, "PCoA", "bray") %>% 
  plot_ordination(ps, ., color='Sequencing_lane', title = "", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Sequencing batch")) +
  theme_bw() + stat_ellipse(level=0.8, lty=2, lwd=0.6) +
  scale_color_manual(values=c('gray50','gray20')) + 
  theme(axis.text = element_text(size = 11))
ordinate(ps, "PCoA", "uni") %>% 
  plot_ordination(ps, ., color='Sequencing_lane', title = "UniFrac", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Sequencing run")) +
  theme_bw() + stat_ellipse(level=0.8, lty=2, lwd=0.6) +
  scale_color_manual(values=c('gray50','gray20')) + 
  theme(axis.text = element_text(size = 11))
ordinate(ps, "PCoA", "wuni") %>% 
  plot_ordination(ps, ., color='Sequencing_lane', title = "Weighted UniFrac", type='samples') +
  geom_point(size=2) +guides(color = guide_legend(title = "Sequencing run")) +
  theme_bw() + stat_ellipse(level=0.8, lty=2, lwd=0.6) +
  scale_color_manual(values=c('gray50','gray20')) + 
  theme(axis.text = element_text(size = 11))
dev.off()


#check total counts
table(ps@sam_data$Species, ps@sam_data$State, ps@sam_data$Year)
table(ps@sam_data$Species)


###-----core again------
core<-data.frame(ps@otu_table) #make a copy of rarefied OTU table
core[core==0]<-NA #turn zeros to NA
core$sumNA<-rowSums(is.na(core)) #add column of sums of NAs per row
summary(core$sumNA) #min is 105, meaning present in (270-105=165) 165/270, or 61% of individuals
core2<-core[which(core$sumNA<135),] #in more than half of the individuals, only 3 OTU
core3<-core[order(core$sumNA),]
core3$prev<-c((270-core3$sumNA)/270) #add column of prevalence (proportion of individuals with this ASV)
core3$sum<-rowSums(core3[,1:270], na.rm=T)
which(core3$sumNA==270) #0
core4<-core3[which(core3$prev>0.3),] #mean prevalence 0.010253

hist(core3$prev, breaks=100, ylim=c(0,100))
summary(core3$prev)
core30<-ps@tax_table[which(rownames(ps@tax_table) %in% rownames(core4)),] #OTUs present in > 30% of individuals
core30<-core30[order(rownames(core30))]
#add prevalence column
prev<-core4$prev[order(rownames(core4))]
core30<-as.data.frame(core30)
core30$prev<-prev
#write.csv(core30, 'core30_2017-2020.csv')


#2020
p20_2<-prune_samples(rownames(sampledata[which(sampledata$Year=='2020'),]), ps)
##make rel abundace plot of p20_2, rarefied 2020 samples only
topp3.11<-as.data.frame(p20_2@tax_table)
sort(table(topp3.11$Phylum), decreasing=T)
topp<-names(sort(table(topp3.11$Phylum), decreasing=T)[1:8]) #top 8 phyla, consistent btwn 2.11 and 3.11

glom3.11<-tax_glom(p20_2, taxrank = "Phylum", NArm=F) #glom by phylum
mglom3.11<-psmelt(glom3.11) 
mglom3.11$Phylum[which(!(mglom3.11$Phylum %in% topp))]<-'other' #change non-common phyla to other 
table(mglom3.11$Phylum)

pdf(file='figs_for_paper/abundance_microbiome_2020.pdf', height=5, width=5)
level_order<-c('BTNW','MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','HOWA','CAWA','COYE','NAWA','BAWW','WEWA','OVEN')
ggplot(mglom3.11, aes(x = factor(sample_Species, level=level_order), y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','turquoise4',
                             'rosybrown3','turquoise2','thistle1',
                             'rosybrown2','wheat4','sandybrown','black')) +
  ggtitle('16s Batch 2') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(legend.position="none")
dev.off()


###-------look at overlap between sequencing runs in OTU ids-------
otu1<-read.csv("~/Documents/Toews_Lab/16s_data/Jan2020/silva/silva_taxonomy.tsv",sep='\t', header=F, skip=2)
otu2<-read.csv("~/Documents/Toews_Lab/16s_2020samples/qiime_workflow/decontam/silva_taxonomy.tsv",sep='\t', header=F, skip=2)

length(unique(otu1$V1)) #7151 OTUs in 2017-2019 dataset
length(unique(otu2$V1)) #12057 OTUs in 2020 dataset

length(otu1$V1[which(otu1$V1 %in% otu2$V1)]) #1164 OTUs from otu1 in otu2
length(otu2$V1[which(otu2$V1 %in% otu1$V1)]) #1164 OTUs from otu2 in otu1
overlap_otus<-(intersect(otu1$V1, otu2$V1)) #confirmed 1164 OTUs in both

venn.diagram(
  x = list(otu1$V1, otu2$V1),
  category.names = c("N OTUs\nLane 1" , "N OTUs\nLane 2"),
  filename = 'otu_venn.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c('coral2','cyan3')
)

##physeq2 overlapping vs unique OTU count histograms
#unique OTUs in 2020 lane
u20<-otu2$V1[which(!(otu2$V1 %in% overlap_otus))]
#unique OTUs in 1st lane
u3<-otu1$V1[which(!(otu1$V1 %in% overlap_otus))]

#get rid of zero sum OTUs
physeq22<-prune_taxa(taxa_sums(physeq2) > 0, physeq2)

summary(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u20)]))
summary(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% overlap_otus)]))
summary(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u3)]))

taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u20)] < 10)

hist(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u20)]), breaks=100, xlim=c(0,10000))
hist(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% overlap_otus)]), add=T, col='blue', breaks=1000)
hist(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u3)]), add=T, col='red', breaks=100)

pdf('density_OTUcounts_zoomed.pdf', height=5, width=5)
plot(density(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u20)])), xlim=c(0,2000),
     xlab='OTU read counts', main='', col='cyan3')
lines(density(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% overlap_otus)])), col='black')
lines(density(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u3)])),col='coral3')
legend('topright', legend=c('2020 unique OTUs', 'Overlapping OTUs', '2017-2019 unique OTUs'),
       col=c('cyan3', 'black', 'coral3'), lty=1, cex=0.9)
dev.off()

pdf('density_OTUcounts.pdf', height=5, width=5)
plot(density(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u20)])),
     xlab='OTU read counts', main='', col='cyan3')
lines(density(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% overlap_otus)])), col='black')
lines(density(taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u3)])),col='coral3')
legend('topright', legend=c('2020 unique OTUs', 'Overlapping OTUs', '2017-2019 unique OTUs'),
       col=c('cyan3', 'black', 'coral3'), lty=1, cex=0.9)
dev.off()

#and rarefied (very similar)
plot(density(taxa_sums(ps@otu_table[which(rownames(ps@otu_table) %in% u20)])),xlim=c(0,1000),
     xlab='OTU read counts', main='Rarefied', col='cyan3')
lines(density(taxa_sums(ps@otu_table[which(rownames(ps@otu_table) %in% overlap_otus)])), col='black')
lines(density(taxa_sums(ps@otu_table[which(rownames(ps@otu_table) %in% u3)])),col='coral3')
legend('topright', legend=c('2020 unique OTUs', 'Overlapping OTUs', '2017-2019 unique OTUs'),
       col=c('cyan3', 'black', 'coral3'), lty=1, cex=0.9)

t1<-taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u20)])
t2<-taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% overlap_otus)])
t3<-taxa_sums(physeq22@otu_table[which(rownames(physeq22@otu_table) %in% u3)])
length(t1[t1<10])
length(t1[t1<50])/length(t1)
length(t2[t2<50])
length(t2[t2<50])/length(t2)
length(t3[t3<50])
length(t3[t3<50])/length(t3)

#find average read counts per indvidiaul
summary(sdt$TotalReads[which(sdt$Sequencing_lane=='1')])
sd(sdt$TotalReads[which(sdt$Sequencing_lane=='1')])
summary(sdt$TotalReads[which(sdt$Sequencing_lane=='2')])
sd(sdt$TotalReads[which(sdt$Sequencing_lane=='2')])

#chao1 and faiths pd
alpha$chao<-estimate_richness(ps, measures=c('Chao1'))[1] #add to alpha dataframe
kruskal.test(alpha$chao$Chao1, sample_data(ps)$Species) #Kruskal-Wallis chi-squared = 13.987, df = 14, p-value = 0.4507
library(btools) #for faiths pd
testpd<-estimate_pd(ps)
alpha$PD<-testpd$PD
kruskal.test(alpha$PD, sample_data(ps)$Species) #Kruskal-Wallis chi-squared =14.979 , df = 14, p-value = 0.3796

#PERMANOVA on jaccard
jac<-distance(ps, method='jaccard', binary=T)
adonis2(jac ~ Species, data=md) #r2=0.06229, p=0.001****
plot(betadisper(jac, md$Species))
permutest(betadisper(jac, md$Species)) #F=4.4809,p=0.001***
adonis2(jac ~ State, data=md) #r2=0.01473, p=0.001****
plot(betadisper(jac, md$State))
permutest(betadisper(jac, md$State)) #F=30.959,p=0.001***
adonis2(jac ~ Year, data=md) #r2=0.09416, p=0.001****
plot(betadisper(jac, md$Year))
permutest(betadisper(jac, md$Year)) #F=151.68,p=0.001***
adonis2(jac ~ Sequencing_lane, data=md) #r2=0.08376, p=0.001****
plot(betadisper(jac, md$Sequencing_lane))
permutest(betadisper(jac, md$Sequencing_lane)) #F=527.04,p=0.001***

adonis2(jac ~ Species+State+Year, data=md, by='margin')
#species r2=0.05263, p=0.004
#state r2=0.00458, p=0.017
#year r2=0.07986, p=0.001
adonis2(jac ~ Species+State+Sequencing_lane, data=md, by='margin')
#species r2=0.05306, p=0.003
#state r2=0.00518, p=0.009
#sequencing land r2=0.07066, p=0.001
table(md$Species, md$State, md$Year, md$Sequencing_lane) #Year and lane do not intersect--has zeros!

adonis2(bray ~ Species+State+Year, data=md, by='margin')
#species r2=0.04705, p=0.192
#state r2=0.00534, p=0.033
#year r2=0.13825, p=0.001
adonis2(bray ~ Species+State+Sequencing_lane, data=md, by='margin')
#species r2=0.04720, p=0.162
#state r2=0.00657, p=0.011
#sequencing lane r2=0.13049, p=0.001

adonis2(uni ~ Species+State+Year, data=md, by='margin')
#species r2=0.05642, p=0.001
#state r2=0.00476, p=0.036
#year r2=0.08698, p=0.001
adonis2(uni ~ Species+State+Sequencing_lane, data=md, by='margin')
#species r2=0.05712, p=0.001
#state r2=0.00551, p=0.012
#sequencing lane r2=0.07931, p=0.001

adonis2(wuni ~ Species+State+Year, data=md, by='margin')
#species r2=0.06124, p=0.180
#state r2=0.00726, p=0.106
#year r2=0.01642, p=0.146
adonis2(wuni ~ Species+State+Sequencing_lane, data=md, by='margin')
#species r2=0.06219, p=0.163
#state r2=0.01286, p=0.020
#sequencing lane r2=0.00572, p=0.194

#export distance matricies for mantel tests 
#full dataset
saveRDS(bray, file='bray_full.rds')
saveRDS(jac, file='jac_full.rds')
saveRDS(uni, file='unifrac_full.rds')
saveRDS(wuni, file='wunifrac_full.rds')

#2020 subset (b2)
jac20<-distance(p20_2, method='jaccard', binary=T)
saveRDS(bray20, file='bray_b2.rds')
saveRDS(jac20, file='jac_b2.rds')
saveRDS(uni20, file='unifrac_b2.rds')
saveRDS(wuni20, file='wunifrac_b2.rds')


#see how many of the 14 species carry each core bacteria
table(ps@sam_data$Species)
head(core30) #core w/prevalence > 30%
amre<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='AMRE')],ps)
amre<-prune_taxa(taxa_sums(amre) > 0, amre) #remove zero sum taxa
amre_c<-as.data.frame(amre@tax_table)
amre_c<-amre_c[which(row.names(amre_c) %in% rownames(core30)),] 
dim(amre_c) #all 39

baww<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='BAWW')],ps)
baww<-prune_taxa(taxa_sums(baww) > 0, baww) #remove zero sum taxa
baww_c<-as.data.frame(baww@tax_table)
baww_c<-baww_c[which(row.names(baww_c) %in% rownames(core30)),] #all 9
dim(baww_c) #all 39

blbw<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='BLBW')],ps)
blbw<-prune_taxa(taxa_sums(blbw) > 0, blbw) #remove zero sum taxa
blbw_c<-as.data.frame(blbw@tax_table)
blbw_c<-blbw_c[which(row.names(blbw_c) %in% rownames(core30)),] #all 9
dim(blbw_c) #all 39

btbw<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='BTBW')],ps)
btbw<-prune_taxa(taxa_sums(btbw) > 0, btbw) #remove zero sum taxa
btbw_c<-as.data.frame(btbw@tax_table)
btbw_c<-btbw_c[which(row.names(btbw_c) %in% rownames(core30)),] #all 9
dim(btbw_c) #all 39

btnw<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='BTNW')],ps)
btnw<-prune_taxa(taxa_sums(btnw) > 0, btnw) #remove zero sum taxa
btnw_c<-as.data.frame(btnw@tax_table)
btnw_c<-btnw_c[which(row.names(btnw_c) %in% rownames(core30)),] #all 9
dim(btnw_c) #all 39

cawa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='CAWA')],ps)
cawa<-prune_taxa(taxa_sums(cawa) > 0, cawa) #remove zero sum taxa
cawa_c<-as.data.frame(cawa@tax_table)
cawa_c<-cawa_c[which(row.names(cawa_c) %in% rownames(core30)),] #all 9
dim(cawa_c) #all 39

coye<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='COYE')],ps)
coye<-prune_taxa(taxa_sums(coye) > 0, coye) #remove zero sum taxa
coye_c<-as.data.frame(coye@tax_table)
coye_c<-coye_c[which(row.names(coye_c) %in% rownames(core30)),] #all 9
dim(coye_c) #all 39

cswa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='CSWA')],ps)
cswa<-prune_taxa(taxa_sums(cswa) > 0, cswa) #remove zero sum taxa
cswa_c<-as.data.frame(cswa@tax_table)
cswa_c<-cswa_c[which(row.names(cswa_c) %in% rownames(core30)),] #all 9
dim(cswa_c) #all 39

howa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='HOWA')],ps)
howa<-prune_taxa(taxa_sums(howa) > 0, howa) #remove zero sum taxa
howa_c<-as.data.frame(howa@tax_table)
howa_c<-howa_c[which(row.names(howa_c) %in% rownames(core30)),] 
dim(howa_c) #all 39 -include all core taxa now, picked up after batch 1

mawa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='MAWA')],ps)
mawa<-prune_taxa(taxa_sums(mawa) > 0, mawa) #remove zero sum taxa
mawa_c<-as.data.frame(mawa@tax_table)
mawa_c<-mawa_c[which(row.names(mawa_c) %in% rownames(core30)),] 
dim(mawa_c) #all 39

mywa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='MYWA')],ps)
mywa<-prune_taxa(taxa_sums(mywa) > 0, mywa) #remove zero sum taxa
mywa_c<-as.data.frame(mywa@tax_table)
mywa_c<-mywa_c[which(row.names(mywa_c) %in% rownames(core30)),] 
dim(mywa_c) #all 39

nawa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='NAWA')],ps)
nawa<-prune_taxa(taxa_sums(nawa) > 0, nawa) #remove zero sum taxa
nawa_c<-as.data.frame(nawa@tax_table)
nawa_c<-nawa_c[which(row.names(nawa_c) %in% rownames(core30)),] 
dim(nawa_c) #only 10!
rownames(core30)[(which(!(rownames(core30) %in% rownames(nawa_c))))]
#"01d458ccec006cf78977428231b8f108" "0278639f012e21f07140177e64a3cd7c"
# "0955475bbf7666f6b1465c81269d37ad" "1e31e490bc9d4e6c44e8c0248240048d"
# "202948c580b39ab15bfe220f01390d98" "2308b53c0d1c7958976087932aacbc13"
# "2f15cf9dea68e788c1c17bb967eb6ff0" "47fbd45f80e40fdc565e816d039d24a9"
# "510282ab2c5fa95af9ef3c68d6f84e6f" "550d54a1386bfb7d9df8076e69906b5a"
# "5a34853b101975f61d8ca4aac3d453fc" "5b3ee3533b213ee4ef18e75127a59e14"
# "65f9076233443330b58ff3a4b8b57826" "6cb83ab369040d99396cf7e24c706354"
# "70697a734288fc8eee4055af01ec84d1" "711634e0bbb04864b5bb168ab4dbfaeb"
# "7d51c3f4b9e05a44231a81eaa0bea39f" "8d716afc0a7075d7fab2fe6bb9a7a544"
# "8e36a3097dcc9c68b0f0b1e70a55452b" "945184b6386c192c0066e0a98a154780"
# "958264b206165f7eb3b34707654d2702" "a7da93a9744a7863292bcf97d47b6a10"
# "b4ea63eb05c7e0ce9ccb755b19d47240" "c9541498e458b40f550b05fe359e951d"
# "cdfb2c9e02e4f387c02256cafdcbeeec" "e6ebb79f2ebe5a66b033c8910bbcf194"
# "ee76384ae44917044fe25a04f0f62792" "fab36fed8422b57daf15bfeaa047aa5b"
# "fb5c2aa7e7a97d1ba29b2c9405cb6a41"

nopa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='NOPA')],ps)
nopa<-prune_taxa(taxa_sums(nopa) > 0, nopa) #remove zero sum taxa
nopa_c<-as.data.frame(nopa@tax_table)
nopa_c<-nopa_c[which(row.names(nopa_c) %in% rownames(core30)),] 
dim(nopa_c) #all 39

oven<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='OVEN')],ps)
oven<-prune_taxa(taxa_sums(oven) > 0, oven) #remove zero sum taxa
oven_c<-as.data.frame(oven@tax_table)
oven_c<-oven_c[which(row.names(oven_c) %in% rownames(core30)),] 
dim(oven_c) #all 39, up from batch 1

wewa<-prune_samples(rownames(ps@sam_data)[which(ps@sam_data$Species=='WEWA')],ps)
wewa<-prune_taxa(taxa_sums(wewa) > 0, wewa) #remove zero sum taxa
wewa_c<-as.data.frame(wewa@tax_table)
wewa_c<-wewa_c[which(row.names(wewa_c) %in% rownames(core30)),] 
dim(wewa_c) #all 39!


#individual relative abundance plot 
level_order2<-data.frame(id=mps2$Sample,species=mps2$sample_Species)
orders<-data.frame(species=level_order, order=1:15)
#merge <- left_join(orders, level_order2, by = "species")
#merge$Sample<-merge$id
#orders[14,2]<-14
mps2$order<-match(mps2$sample_Species, orders$species)
mps2$order2<-paste(mps2$order, mps2$Sample, sep='_')
mps2$new<-mps2$order
mps2$new <- car::recode(mps2$new, '10=91;11=92;12=93;13=94;14=95;15=96')
table(mps2$new) 
table(mps2$sample_Species)
mps2$order3<-paste(mps2$new, mps2$Sample, sep='_')

pdf(file='~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/figs/ind_abund_micro.pdf', height=4, width=12)
ggplot(mps2, aes(x = order3, y = Abundance, fill = Phylum2)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','turquoise4',
                             'rosybrown3','thistle1','violetred4',
                             'rosybrown2','wheat4','sandybrown','black')) +
  ggtitle('16s') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(axis.text = element_text(size = 11)) +
  theme(legend.position="none") + geom_vline(xintercept=26.5 ) + geom_vline(xintercept=46.5 ) +
  geom_vline(xintercept=64.5 ) + geom_vline(xintercept=84.5 ) + geom_vline(xintercept=101.5 ) + 
  geom_vline(xintercept=117.5 ) + geom_vline(xintercept=145.5 ) + geom_vline(xintercept=168.5 ) +
  geom_vline(xintercept=179.5 ) + geom_vline(xintercept=196.5 ) + geom_vline(xintercept=210.5 ) +
  geom_vline(xintercept=221.5 ) + geom_vline(xintercept=243.5 ) + geom_vline(xintercept=248.5 ) +
  geom_vline(xintercept=270.5) + theme(axis.text.x = element_blank())
dev.off()

table(md$Species)

#how many not identified beyond bacteria?
tt<-data.frame(ps@tax_table)
table(tt$Phylum, useNA = 'always')
