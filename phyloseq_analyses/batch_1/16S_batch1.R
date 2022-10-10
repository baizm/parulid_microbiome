setwd("~/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/")
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
load('phyloseq4.RData')

#read in post-decontam filtered data from qiime2
########-------create physeq object---------------------########
####------read in OTU table, has to be a matrix---------------------
otu_table<-read.csv('feature-table.tsv', sep='\t', header=T, skip=1)
#rearrange columns so negatives are at end to match sampledata (don't know if necessary)
otu_table<-otu_table[,c(1,4:253,2:3)]
#convert to matrix
otumat<-as.matrix(otu_table[,2:253])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
#correct negative names so they match other tables
colnames(otumat)[251:252]<-c('7-neg-extr', '8-neg-extr')
###-------read in taxonomy table-----------------
tax_table<-read.csv('../../silva_taxonomy.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#it wasn't filtered, so i'll have to delete mt,cp,un,euk,decontam...
tax_table2<-tax_table2[which(tax_table2$Feature.ID %in% otu_table$X.OTU.ID),]
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('../../../metadata2_16s_ny.tsv', sep='\t', header=T,)
rownames(sampledata)<-sampledata$id
x<-sampledata$id[(which(!(sampledata$id %in% colnames(otumat))))]
sampledata<-sampledata[which(!(sampledata$id %in% x)),] #get rid of 3 samples filtered by qiime (2 neg + 1 low)
sampledata<-sampledata[,2:29]
###-------read in nwk tree with ape package-----------
tree<-read.tree('tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
physeq = phyloseq(OTU, TAX, SAM, tree)
#sums of OTUs across samples
taxa_sums(physeq)
#alpha - tells me I don't have singletons...is this suspicious?
plot_richness(physeq, x="Species", color="Species",measures="Shannon") + theme(legend.position="none") 


#uni<-distance(physeq, method='unifrac') ###issues STILLLLLL
##use physeq to calculate CORE microbiome, after filtering out <4k depth samples, below

#write file of pruned warblers for qiime filtering
#prune species with < 4 indiv 
l4<-c('BLPW','BWWA','GWWA','LOWA','NOWA','PIWA','YEWA') #species with less than 5 indiv
#new pruned dataset removing poorly sampled species
physeq2<-prune_samples(!(grepl(paste(l4,collapse='|'),sample_data(physeq)$Species)), physeq)
write.csv(rownames(physeq2@sam_data), './qiime_diversity/samples_to_keep.txt',quote=F,row.names=F)



####---------find samples to remove due to low read counts & negatives-------------
#make data frame with total read counts after filtering using sample_sums
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
summary(sdt$TotalReads)
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
#look at new rarecurve
#pdf(file='rarecurve.pdf', height=4, width=4)
rarecurve(t(otu_table(physeq)), step=50, lwd=0.8, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,9000))
#dev.off()
#data frame of samples with reads below threshold
low<-sdt[which(sdt$TotalReads<4000),]
table(low$Species) 
table(physeq@sam_data$Species) #for comparison, it looks like mostly removes from species that already have a lot of individuals (except GWWA)
sdt[which(sdt$Species=='GWWA'),] #but the two have ridiculously low read counts anyway
table(low$Type) #includes neg controls and 81 pos samples

#make vector of samples to keep
keep_samples<-sdt$SampleID[which(!(sdt$SampleID %in% low$SampleID))]
#make new phyloseq object to exclude low read samples
#ps2 is non-rarefied dataset! 169 individuals
ps2<-prune_samples(keep_samples, physeq)
#plot original vs pruned samples by species
par(mfrow=c(2,1), mar=c(4,4,2,1))
barplot(table(physeq@sam_data$Species), las=2, ylab='N', main='Original')
abline(h=4, col='red')
barplot(table(ps2@sam_data$Species), las=2, ylab='N', main='Pruned')
abline(h=4, col='red')
min(sample_sums(ps2)) #4020
max(sample_sums(ps2)) #155961


###--------identify the core microbiome, general warbler GM -----------
core<-data.frame(ps2@otu_table) #make a copy of non-rarefied OTU table
core[core==0]<-NA #turn zeros to NA
core$sumNA<-rowSums(is.na(core)) #add column of sums of NAs per row
summary(core$sumNA) #min is 66, meaning present in (169-66=103) 103/169, or 61% of individuals
core2<-core[which(core$sumNA<84),] #in more than half of the individuals, only 1 OTU: 63afe8e6aac58bf0d670a82ca5bc574c
core3<-core[order(core$sumNA),]
core3$prev<-c((169-core3$sumNA)/169) #add column of prevalence (proportion of individuals with this ASV)
core3$sum<-rowSums(core3[,1:169], na.rm=T)
which(core3$sumNA==169)


save.image('phyloseq4.RData')


##do in qiime2----------------
###create rarefied data set=ps3
#find lowest N reads for threshold
sdt2 = data.table(as(sample_data(ps2), "data.frame"),TotalReads = sample_sums(ps2), keep.rownames = TRUE)
min(sdt2$TotalReads) #4020
ps3<-rarefy_even_depth(ps2,sample.size=4020,rngseed = 999, replace = F, trimOTUs = T, verbose = T)

###--------identify the core microbiome, general warbler GM -----------
core<-data.frame(ps2@otu_table) #make a copy of non-rarefied OTU table
core[core==0]<-NA #turn zeros to NA
core$sumNA<-rowSums(is.na(core)) #add column of sums of NAs per row
summary(core$sumNA) #min is 66, meaning present in (169-66=103) 103/169, or 61% of individuals
core2<-core[which(core$sumNA<84),] #in more than half of the individuals, only 1 OTU: 63afe8e6aac58bf0d670a82ca5bc574c
core3<-core[order(core$sumNA),]
core3$prev<-c((169-core3$sumNA)/169) #add column of prevalence (proportion of individuals with this ASV)
core3$sum<-rowSums(core3[,1:169], na.rm=T)
which(core3$sumNA==169)
zo1<-which(core3$sum==0) #same as above. indexes of rows, remove all these OTUs
ro1<-rownames(core3[zo1,]) #otus to remove from ps2.1 bc prevalence is zero

hist(core3$prev,breaks=50, xlab="prevalence across individuals",main='')
hist(core3$prev,breaks=50, xlim=c(0.1,0.62), ylim=c(0,30),
     xlab="prevalence across individuals",main='') #zoomed in

#prune species with < 4 indiv 
l4<-c('BLPW','BWWA','GWWA','LOWA','NOWA','PIWA','YEWA') #species with less than 4 indiv
#new pruned datasets removing poorly sampled species
ps2.1<-prune_samples(!(grepl(paste(l4,collapse='|'),sample_data(ps2)$Species)), ps2)
ps3.1<-prune_samples(!(grepl(paste(l4,collapse='|'),sample_data(ps3)$Species)), ps3)

###prune OTUs with zero sums
#non-rarefied
ro1.2<-which(!(taxa_names(ps2.1) %in% ro1))
ro1.22<-taxa_names(ps2.1)[c(ro1.2)]
#make sure none of the 89 ASVs are in v2
ps2.11<-prune_taxa(ro1.22, ps2.1) #ps2.11 is pruned warbler spp and pruned OTUs

#rarefied
core_r<-data.frame(ps3.1@otu_table) #make a copy of rarefied OTU table
core_r[core_r==0]<-NA #turn zeros to NA
core_r$sumNA<-rowSums(is.na(core_r)) #add column of sums of NAs per row
summary(core_r$sumNA) #min is 65
core2_r<-core_r[which(core_r$sumNA<80),] #in more than half of the individuals
core3_r<-core_r[order(core_r$sumNA),]
core3_r$prev<-c((161-core3_r$sumNA)/161) #add column of prevalence (proportion of individuals with this ASV)
core3_r$sum<-rowSums(core3_r[,1:161], na.rm=T)
which(core3_r$sumNA==161)
zo2<-which(core3_r$sum==0) #same as above. indexes of rows, remove all these OTUs
ro2<-rownames(core3_r[zo2,]) #otus to remove from ps3.1 bc prevalence is zero
ro2.2<-which(!(taxa_names(ps3.1) %in% ro2))
ro2.22<-taxa_names(ps3.1)[c(ro2.2)]
#make sure none of the 89 ASVs are in v2
ps3.11<-prune_taxa(ro2.22, ps3.1) #ps3.11 is pruned warbler spp and pruned OTUs


#######alpha diversity#######
#alpha diversity rarefied, with poorly sampled species removed
ar<-plot_richness(ps3.11, x='Species', color='Species',measures=c('Observed','Shannon')) + 
  geom_boxplot() +
  theme(legend.position="none") 
#alpha diversity non-rarefied
an<-plot_richness(ps2.11, x='Species', color='Species',measures=c('Observed','Shannon')) + 
  geom_boxplot() +
  theme(legend.position="none") 
#pdf(file='alpha_comparison.pdf', height=6.5, width=5.5)
grid.arrange(an+ggtitle("non-rarefied"),ar+ggtitle('rarefied'),ncol=1) #alpha div for rare-vs-non-rarefied look very similar!
#dev.off()

#calculate alpha values
alpha_r<-estimate_richness(ps3.11, measures=c('Shannon','Observed'))
kruskal.test(alpha_r$Observed, sample_data(ps3.11)$Species) #Kruskal-Wallis chi-squared = 20.115, df = 13, p-value = 0.09237



#plot distinguishing only top 10 phyla
topp2.11<-as.data.frame(ps2.11@tax_table)
sort(table(topp2.11$Phylum), decreasing=T)
topp3.11<-as.data.frame(ps3.11@tax_table)
sort(table(topp3.11$Phylum), decreasing=T)
topp<-names(sort(table(topp3.11$Phylum), decreasing=T)[1:8]) #top 8 phyla, consistent btwn 2.11 and 3.11

#taxa bar plot for non-rarefied data, top phyla
glom2.11<-tax_glom(ps2.11, taxrank = "Phylum", NArm=F) #glom by phylum
mglom2.11<-psmelt(glom2.11)
mglom2.11$Phylum[which(!(mglom2.11$Phylum %in% topp))]<-'other' #change non-common phyla to other 
table(mglom2.11$Phylum)
#pdf(file='taxaPlot_non-rarefied.pdf', height=6, width=8)
ggplot(mglom2.11, aes(x = sample_Species, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('cornsilk3','darkblue',
                             'cornflowerblue','coral3','darkgreen',
                             'darkolivegreen3','black','darkgoldenrod','gold2')) +
  ggtitle('Non-rarefied dataset') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
#dev.off()
#and for rarefied dataset
glom3.11<-tax_glom(ps3.11, taxrank = "Phylum", NArm=F) #glom by phylum
mglom3.11<-psmelt(glom3.11) #39 OTUs x 161 individuals=6279 rows
mglom3.11$Phylum[which(!(mglom3.11$Phylum %in% topp))]<-'other' #change non-common phyla to other 
table(mglom3.11$Phylum)
#pdf(file='taxaPlot_rarefied.pdf', height=6, width=8)
ggplot(mglom3.11, aes(x = sample_Species, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('cornsilk3','darkblue',
                             'cornflowerblue','coral3','darkgreen',
                             'darkolivegreen3','black','darkgoldenrod','gold2')) +
  ggtitle('Rarefied dataset') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
#dev.off()

####-----------beta diversity-----------------------------
#test for effect of year and state to see if I can pool
#calculate distance matrix non-rare
bray2.11<-distance(ps2.11, method='bray')
#and unifrac
uni2.11<-distance(ps2.11, method='unifrac')
#and weighted unifrac
wuni2.11<-distance(ps2.11, method='wunifrac')

##is error bc of multibranching?? https://github.com/joey711/phyloseq/issues/936
edges=phy_tree(ps3.11)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
#[1] 4970
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
#[1] 6
mycounts[mycounts !=2] # How many nodes each of the above has
#no...can't figure it out. Export back to qiime?
#to do once in qiime
#rarefy
#remove contaminant ASVs
#create subset of pruned warbler spp.
#diversity stats




#make dataframe from sample data
md2.1<-data.frame(sample_data(ps2.1), row.names=rownames(sample_data(ps2.1)))
md2.1$Year<-as.factor(md2.1$Year) #make non-continuous
ps2.1@sam_data$Year<-as.factor(ps2.1@sam_data$Year) #make non-continuous
md2.1$doy<-yday(as.Date(md2.1$Date, '%m/%d/%y')) #make day of year column
#adonis tests, non-rare
a1<-adonis(bray2.1 ~ Year*Species, data=md2.1) #non-significant interaction, can combine years?
a2<-adonis(bray2.1 ~ State*Species, data=md2.1) #interaction significant, cannot combine NY-PA?
dis1<-permutest(betadisper(bray2.1, md2.1$Year)) #non-significant, adonis result not due to differences in group dispersions?
dis2<-permutest(betadisper(bray2.1, md2.1$State)) #non-significant, adonis result not due to differences in group dispersions?
a3<-adonis(uni2.1 ~ Year*Species, data=md2.1) #non-significant interaction, can combine years?
a4<-adonis(uni2.1 ~ State*Species, data=md2.1) #non-significant, can combine NY-PA?
dis3<-permutest(betadisper(uni2.1, md2.1$Year)) #significant, adonis result could be due to differences in group dispersions? But year not sig....so doesn't matter?
dis4<-permutest(betadisper(uni2.1, md2.1$State)) #non-significant, adonis result not due to differences in group dispersions?
a5<-adonis(wuni2.1 ~ Year*Species, data=md2.1) #non-significant interaction, can combine years?
a6<-adonis(wuni2.1 ~ State*Species, data=md2.1) #not-significant interaction, can combine NY-PA?
dis5<-permutest(betadisper(wuni2.1, md2.1$Year)) #non-significant, adonis result not due to differences in group dispersions? 
dis6<-permutest(betadisper(wuni2.1, md2.1$State)) #non-significant, adonis result not due to differences in group dispersions?
#write to file output of test for year, state, non-rare
beta_out1<-capture.output(print(a1),print(a2),print(dis1),print(dis2),
                          print(a3),print(a4),print(dis3),print(dis4),
                          print(a5),print(a6),print(dis5),print(dis6))
writeLines(beta_out1, con = file("result_NR_YearState.txt"))

##calculate distance matrix rarefied
bray3.11<-distance(ps3.11, method='bray')
#and unifrac
uni3.11<-distance(ps3.11, method='unifrac')
#and weighted unifrac
wuni3.11<-distance(ps3.11, method='wunifrac')
#make dataframe from sample data
md3.11<-data.frame(sample_data(ps3.11), row.names=rownames(sample_data(ps3.11)))
md3.11$Year<-as.factor(md3.11$Year) #make non-continuous
ps3.11@sam_data$Year<-as.factor(ps3.11@sam_data$Year) #make non-continuous
md3.11$doy<-yday(as.Date(md3.11$Date, '%m/%d/%y')) #make day of year column
plot(alpha_r$Observed~md3.11$doy, col=factor(md3.11$Species))
pdf(file='~/Documents/Academic stuff/faculty_interview_stuff/BoiseState/temporal.pdf', height=4, width=5)
plot(alpha_r$Observed~md3.11$doy, xlab='Day', ylab='Diversity', lwd=2)
dev.off()
fit <- lm(alpha_r$Observed~poly(md3.11$doy,2,raw=TRUE))
summary(fit)

md3.11$temporal[md3.11$doy < 160]<-'early'
md3.11$temporal[md3.11$doy > 170]<-'late'
md3.11$temporal[md3.11$doy <= 170 & md3.11$doy >= 160]<-'mid'
md3.11$observed<-alpha_r$Observed
boxplot(md3.11$observed~md3.11$temporal)
summary(aov(md3.11$observed~md3.11$temporal))
#try coding mid point based on each species max observed date

#adonis tests, rare
a7<-adonis(bray3.11 ~ Year*Species, data=md3.11) #non-significant interaction, can combine years?
a8<-adonis(bray3.11 ~ State*Species, data=md3.11) #interaction significant, cannot combine NY-PA?
dis7<-permutest(betadisper(bray3.11, md3.11$Year)) #non-significant, adonis result not due to differences in group dispersions?
dis8<-permutest(betadisper(bray3.11, md3.11$State)) #non-significant, adonis result not due to differences in group dispersions?
a9<-adonis(uni3.11 ~ Year*Species, data=md3.11) #non-significant interaction, can combine years?
a10<-adonis(uni3.11 ~ State*Species, data=md3.11) #non-significant, can combine NY-PA?
dis9<-permutest(betadisper(uni3.11, md3.11$Year)) #non-significant
dis10<-permutest(betadisper(uni3.11, md3.11$State)) #non-significant, adonis result not due to differences in group dispersions?
a11<-adonis(wuni3.11 ~ Year*Species, data=md3.11) #non-significant interaction, can combine years?
a12<-adonis(wuni3.11 ~ State*Species, data=md3.11) #not-significant interaction, can combine NY-PA?
dis11<-permutest(betadisper(wuni3.11, md3.11$Year)) #non-significant, adonis result not due to differences in group dispersions? 
dis12<-permutest(betadisper(wuni3.11, md3.11$State)) #non-significant, adonis result not due to differences in group dispersions?
#write to file output of test for year, state, non-rare
beta_out2<-capture.output(print(a7),print(a8),print(dis7),print(dis8),
                          print(a9),print(a10),print(dis9),print(dis10),
                          print(a11),print(a12),print(dis11),print(dis12))
writeLines(beta_out2, con = file("result_R_YearState.txt"))
 
#use unname(matrix) or adonis2 to avoid error!
adonis(bray3.11 ~ Species, data=md3.11) #r2=0.09282, p=0.022***
plot(betadisper(bray3.11, md3.11$Species), hull=F, ellipse=T, col=pal)
permutest(betadisper(bray3.11, md3.11$Species)) #ns F=0.703, p=0.762
adonis(bray3.11 ~ State, data=md3.11) #r2=0.01881, p=0.001***
permutest(betadisper(bray3.11, md3.11$State)) #ns F=0.7587, p=0.397
adonis(bray3.11 ~ Year, data=md3.11) #r2=0.01545, p=0.081
permutest(betadisper(bray3.11, md3.11$Year)) #ns p=0.17

adonis(uni3.11 ~ Species, data=md3.11) #r2=0.10604, p=0.001***
plot(betadisper(uni3.11, md3.11$Species))
plot(betadisper(uni3.11, md3.11$Species), hull=F, ellipse=T)
permutest(betadisper(uni3.11, md3.11$Species)) # F=3.0492, p=0.002***
adonis(uni3.11 ~ State, data=md3.11) #r2=0.0115, p=0.005***
plot(betadisper(uni3.11, md3.11$State))
permutest(betadisper(uni3.11, md3.11$State)) #ns F=0.0019, p=0.971
adonis(uni3.11 ~ Year, data=md3.11) #r2=0.01618, p=0.027***
permutest(betadisper(uni3.11, md3.11$Year)) #ns F=2.2121, p=0.099

adonis(wuni3.11 ~ Species, data=md3.11) #r2=0.08393, p=0.414
permutest(betadisper(wuni3.11, md3.11$Species)) #p=0.004***
adonis(wuni3.11 ~ State, data=md3.11) #r2=0.01744, p=0.048***
permutest(betadisper(wuni3.11, md3.11$State)) #ns F=0.1831, p=0.655
adonis(wuni3.11 ~ Year, data=md3.11) #r2=0.01555, p=0.292
permutest(betadisper(wuni3.11, md3.11$Year)) #ns p=0.918

#make supp fig of betadispers
pdf(file='betadisper_1_2017-2019.pdf', height=5, width=5)
plot(betadisper(bray3.11, md3.11$Species), hull=F, ellipse=T, col=pal,main='Bray-Curtis',
     seg.lwd = 0.5,sub = '',label.cex=0.5)
dev.off()
col.fill.rect <- addAlpha(col2rgb("white"), alpha=0.9)
#see end for custom function to change label color
pdf(file='betadisper_1_2017-2019_v2.pdf', height=5, width=5)
myplotbetadisper(betadisper(bray3.11, md3.11$Species), ellipse = TRUE, hull = FALSE, 
                 fillrect=col.fill.rect, seg.lwd = 0.5,
                 alphaPoints=0.7, labPoints=NULL, sub='',
                 main= "Bray-Curtis", coltextrect=pal, label.cex=0.6)
dev.off()
pdf(file='betadisper_2_2017-2019.pdf', height=5, width=5)
plot(betadisper(uni3.11, md3.11$Species), hull=F, ellipse=T, col=pal, main='UniFrac',
     seg.lwd = 0.5,sub = '',label.cex=0.5)
dev.off()
pdf(file='betadisper_2_2017-2019_v2.pdf', height=5, width=5)
myplotbetadisper(betadisper(uni3.11, md3.11$Species), ellipse = TRUE, hull = FALSE, 
                 fillrect=col.fill.rect, seg.lwd = 0.5,
                 alphaPoints=0.7, labPoints=NULL, sub='',
                 main= "UniFrac", coltextrect=pal, label.cex=0.6)
dev.off()


#check total counts
table(ps3.11@sam_data$Species, ps3.11@sam_data$State, ps3.11@sam_data$Year)
table(ps3.11@sam_data$Species)
#make genus column for md's
species<-unique(md2.1$Species) #same for md3.1
genus<-c('Setophaga','Leiothlypis','Setophaga','Setophaga','Setophaga','Setophaga','Setophaga',
         'Setophaga','Setophaga','Cardellina','Geothlypis','Seiurus','Mniotilta','Setophaga')
bt<-data.frame(Species=species, Genus=genus)
md2.1<-join(md2.1,bt,by='Species')
md3.1<-join(md3.1,bt,by='Species')
#check if I can copy rownames, and genus column to ps2
identical(md2.1$Band,ps2.1@sam_data$Band) #TRUE
identical(md3.1$Band,ps3.1@sam_data$Band) #TRUE
rownames(md2.1)<-rownames(ps2.1@sam_data)
rownames(md3.1)<-rownames(ps3.1@sam_data)
ps2.1@sam_data$Genus<-md2.1$genus
ps3.1@sam_data$Genus<-md3.1$genus



#look at ordination plot
library(randomcoloR)
#palette<-distinctColorPalette(14)
pal<-palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
          "red2", "#75DED6", "#75E563", "black", "#B44AE1"))
ordinate(ps3.11, "PCoA", "bray") %>% 
  plot_ordination(ps3.11, ., color='Species', shape = "State", title = "Bray-Curtis (rarefied)", type='samples') +
  scale_color_manual(values=c('antiquewhite4', 'aquamarine4', 'azure4','black','blue4',
                              'brown1','burlywood','chartreuse4','cornflowerblue','darkgoldenrod',
                              'darkgray','darkviolet','khaki','lightpink4')) +
  geom_point(size=2) 

#with color palette
#check how many females:
table(ps3.11@sam_data$Sex) #M=148, F=3, sex?=10
ord_b<-ordinate(ps3.11, "PCoA", "bray")

pdf(file='bray_pcoa2017-2019.pdf', height=5, width=5)
plot_ordination(ps3.11, ord_b, color='Species', shape = "State", title = "", type='samples') +
  scale_color_manual(values=pal) + geom_point(size=3) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(shape="Locality", colour="Species")
dev.off()
#ordinate(ps3.11, "PCoA", "bray") %>% 
#  plot_ordination(ps3.11, ., color='Species', shape = "State", title = "Bray-Curtis (rarefied)", type='samples') +
#  scale_color_manual(values=palette) + geom_point(size=3) +
#  theme_bw()
ord_u<-ordinate(ps3.11, "PCoA", "unifrac")
plot_ordination(ps3.11, ord_u, color='Species', shape = "State", title = "Unifrac (rarefied)", type='samples') +
  scale_color_manual(values=pal) + geom_point(size=3) +
  theme_bw()
#ordinate(ps3.11, "PCoA", "unifrac") %>% 
#  plot_ordination(ps3.11, ., color='Species', shape = "State", title = "Unifrac (rarefied)", type='samples') +
#  scale_color_manual(values=palette) + geom_point(size=3) +
#  theme_bw()
ordinate(ps3.11, "PCoA", "wunifrac") %>% 
  plot_ordination(ps3.11, ., color='Species', shape = "State", title = "Weighted Unifrac (rarefied)", type='samples') +
  scale_color_manual(values=palette) + geom_point(size=3) +
  theme_bw()


ordinate(ps2.1, "PCoA", "bray") %>% 
  plot_ordination(ps2.1, ., color='Species', shape = "State", title = "Bray-Curtis (non-rarefied)", type='samples') +
  scale_color_manual(values=c('antiquewhite4', 'aquamarine4', 'azure4','black','blue4',
                              'brown1','burlywood','chartreuse4','cornflowerblue','darkgoldenrod',
                              'darkgray','darkviolet','khaki','lightpink4')) +
  geom_point(size=2)
#+ stat_ellipse()
#+scale_shape_manual(values=c(1,2)) + geom_point(aes(size=Year))

##some test imply effect of state, so subset dataset and re-test with even sampling



#playing around with merging to get average OTU abundance by species, see https://joey711.github.io/phyloseq/merge.html
test<-merge_samples(ps2.1, "Species")

####-----take a closer look at effect of year by subsetting----------------------------------
###look at species sampled well across each year (only in NY, at least 2 each year), rarefied
table(md3.11$Species,md3.11$Year,md3.11$State) #AMRE,BTBW,BTNW,CAWA,MYWA,NAWA,NOPA,OVEN
yrk<-c('AMRE','BTBW','BTNW','CAWA','MYWA','NAWA','NOPA','OVEN') #spp with at least indv in each year
#yrk<-c('MYWA') #look at one at a time
yr<-md3.11[which(md3.11$Species %in% yrk & md3.11$State=="NY"),]
psyr<-prune_samples(rownames(yr), ps3.11)
#looks like years overlap generally
ordinate(psyr, "PCoA", "bray") %>% 
  plot_ordination(psyr, ., color='Species', shape = "Year", title = "Year Bray-Curtis", type='samples') +
  scale_color_manual(values=palette)+
  geom_point(size=3)
yr1<-md3.11[which(md3.11$Species=="AMRE" & md3.11$State=="NY"),]
psyr1<-prune_samples(rownames(yr1), ps3.11)
yr2<-md3.11[which(md3.11$Species=="BTBW" & md3.11$State=="NY"),]
psyr2<-prune_samples(rownames(yr2), ps3.11)
yr3<-md3.11[which(md3.11$Species=="BTNW" & md3.11$State=="NY"),]
psyr3<-prune_samples(rownames(yr3), ps3.11)
yr4<-md3.11[which(md3.11$Species=="CAWA" & md3.11$State=="NY"),]
psyr4<-prune_samples(rownames(yr4), ps3.11)
yr5<-md3.11[which(md3.11$Species=="MYWA" & md3.11$State=="NY"),]
psyr5<-prune_samples(rownames(yr5), ps3.11)
yr6<-md3.11[which(md3.11$Species=="NAWA" & md3.11$State=="NY"),]
psyr6<-prune_samples(rownames(yr6), ps3.11)
yr7<-md3.11[which(md3.11$Species=="NOPA" & md3.11$State=="NY"),]
psyr7<-prune_samples(rownames(yr7), ps3.11)
yr8<-md3.11[which(md3.11$Species=="OVEN" & md3.11$State=="NY"),]
psyr8<-prune_samples(rownames(yr8), ps3.11)

#plot separate species
y1<-ordinate(psyr1, "PCoA", "bray") %>% 
  plot_ordination(psyr1, ., color='Year', shape = "Species", title = "AMRE", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y2<-ordinate(psyr2, "PCoA", "bray") %>% 
  plot_ordination(psyr2, ., color='Year', shape = "Species", title = "BTBW", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y3<-ordinate(psyr3, "PCoA", "bray") %>% 
  plot_ordination(psyr3, ., color='Year', shape = "Species", title = "BTNW", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y4<-ordinate(psyr4, "PCoA", "bray") %>% 
  plot_ordination(psyr4, ., color='Year', shape = "Species", title = "CAWA", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y5<-ordinate(psyr5, "PCoA", "bray") %>% 
  plot_ordination(psyr5, ., color='Year', shape = "Species", title = "MYWA", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y6<-ordinate(psyr6, "PCoA", "bray") %>% 
  plot_ordination(psyr6, ., color='Year', shape = "Species", title = "NAWA", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y7<-ordinate(psyr7, "PCoA", "bray") %>% 
  plot_ordination(psyr7, ., color='Year', shape = "Species", title = "NOPA", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
y8<-ordinate(psyr8, "PCoA", "bray") %>% 
  plot_ordination(psyr8, ., color='Year', shape = "Species", title = "OVEN", type='samples') +
  scale_color_manual(values=c('deeppink', 'cornflowerblue','black')) + geom_point(size=4)
pdf(file='year_sppNY_rarefied.pdf', height=8, width=8)
grid.arrange(y1,y2,y3,y4,y5,y6,y7,y8,ncol=2)
dev.off()


#abundance plot for AOS
#pdf(file='~/Documents/Conferences:Presentations/AOS-SCO2021/abundance_microbiome.pdf', height=5, width=5)
pdf(file='abundance_microbiome_2017-2019.pdf', height=5, width=5.2)
level_order<-c('BTNW','MYWA','CSWA','BLBW','MAWA','NOPA','BTBW','AMRE','HOWA','CAWA','COYE','NAWA','BAWW','WEWA','OVEN')
ggplot(mglom3.11, aes(x = factor(sample_Species, level=level_order), y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','turquoise4',
                             'rosybrown3','thistle1','violetred4',
                             'rosybrown2','wheat4','sandybrown','black')) +
  ggtitle('16s Batch 1') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(legend.position="none")
dev.off()

#phylogenetic div
library(metagMisc)
library(PhyloMeasures)
micro_pd<-phyloseq_phylo_div(ps3.11, measures='PD')
alpha_r$PD<-micro_pd$PD
kruskal.test(alpha_r$Observed, sample_data(ps3.11)$Species) #Kruskal-Wallis chi-squared = 20.115, df = 13, p-value = 0.09237
kruskal.test(alpha_r$Shannon, sample_data(ps3.11)$Species) #Kruskal-Wallis chi-squared = 16.354, df = 13, p-value = 0.2305
kruskal.test(alpha_r$PD, sample_data(ps3.11)$Species) #Kruskal-Wallis chi-squared = 18.764, df = 13, p-value = 0.1306

###-----core again maybe?_------
core<-data.frame(ps3.11@otu_table) #make a copy of rarefied OTU table
core[core==0]<-NA #turn zeros to NA
core$sumNA<-rowSums(is.na(core)) #add column of sums of NAs per row
summary(core$sumNA) #min is 65, meaning present in (161-65=96) 96/161, or 60% of individuals
core2<-core[which(core$sumNA<81),] #in more than half of the individuals, only 1 OTU: 63afe8e6aac58bf0d670a82ca5bc574c
core3<-core[order(core$sumNA),]
core3$prev<-c((161-core3$sumNA)/161) #add column of prevalence (proportion of individuals with this ASV)
core3$sum<-rowSums(core3[,1:161], na.rm=T)
which(core3$sumNA==161) #0
core4<-core3[which(core3$prev>0.3),] #mean prevalence core3 is 0.012710
summary(core4$prev)

#add prevalence column to tax table
prevtab<-as.data.frame(ps3.11@tax_table)
prevtab$prev <- core3$prev[match(rownames(prevtab), rownames(core3))]
#find clostridia
clost<-prevtab[which(prevtab$Class=='c__Clostridia'),]
clost<-clost[order(clost$prev, decreasing = T),] #all are <10% prevalence


hist(core3$prev, breaks=100, ylim=c(0,100))
summary(core3$prev)
core30<-ps3.11@tax_table[which(rownames(ps3.11@tax_table) %in% rownames(core4)),] #OTUs present in > 30% of individuals
#write.csv(core30, 'core30_2017-2019.csv')


####-------------take a closer look at state by subsetting------------------

save.image('phyloseq4.RData')







####=----------functions-----------------
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#####myplotbetadisp.r
###from https://pastebin.com/qhjvZmhB
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
    myordilabel(x, display = "centroids", choices = axes, cex = label.cex, fill=fillrect, col=coltextrect)
  }
  localTitle(main = main, xlab = xlab, ylab = ylab, sub = sub, 
             ...)
  localAxis(1, ...)
  localAxis(2, ...)
  localBox(...)
  class(g) <- "ordiplot"
  invisible(g)
}


myordilabel <- function (x, display, labels, choices = c(1, 2), priority, select, 
                         cex = 0.8, fill = "white", border = NULL, col = NULL, xpd = TRUE, 
                         ...) 
{
  if (missing(display)) 
    display <- "sites"
  x <- scores(x, choices = choices, display = display, ...)
  if (missing(labels)) 
    labels <- rownames(x)
  if (!missing(select)) {
    x <- .checkSelect(select, x)
    labels <- .checkSelect(select, labels)
  }
  if (!missing(priority)) {
    if (!missing(select)) 
      priority <- priority[select]
    ord <- order(priority)
    x <- x[ord, ]
    labels <- labels[ord]
  }
  else {
    ord <- seq_along(labels)
  }
  em <- strwidth("m", cex = cex, ...)
  ex <- strheight("x", cex = cex, ...)
  w <- (strwidth(labels, cex = cex, ...) + em/1.5)/2
  h <- (strheight(labels, cex = cex, ...) + ex/1.5)/2
  if (is.null(col)) 
    if (!is.null(border)) 
      col <- border
  else col <- par("fg")
  col <- rep(col, length = nrow(x))[ord]
  if (!is.null(border)) 
    border <- rep(border, length = nrow(x))[ord]
  fill <- rep(fill, length = nrow(x))[ord]
  for (i in 1:nrow(x)) {
    vegan:::ordiArgAbsorber(x[i, 1] + c(-1, 1, 1, -1) * w[i], x[i, 
                                                                2] + c(-1, -1, 1, 1) * h[i], col = fill[i], border = border[i], 
                            xpd = xpd, FUN = polygon, ...)
    vegan:::ordiArgAbsorber(x[i, 1], x[i, 2], labels = labels[i], 
                            cex = cex, col = col[i], xpd = xpd, FUN = text, ...)
  }
  invisible(x)
}

addAlpha <- function(col, alpha) {
  alpha <- round(alpha*255)
  rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
}


####for reviews May 2022
pal<-palette(c("#6D7BD0", "blue1", "#D8C278", "#D5E356", "#D45691", "#D988DF", "#D4E1CE", "#79ACCD", "#C1958E",
               "red2", "#75DED6", "#75E563", "black", "#B44AE1"))
ordinate(ps3.11, "NMDS", "bray") %>% 
  plot_ordination(ps3.11, ., color='Species', shape = "State", title = "Bray-Curtis (rarefied)", type='samples') +
  scale_color_manual(values=pal) +
  geom_point(size=3) 

plot_ordination(ps3.11, ord_b, color='Species', shape = "State", title = "", type='samples',
                axes=c(1,3)) +
  scale_color_manual(values=pal) + geom_point(size=3) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(shape="Locality", colour="Species")

#find mean asvs per individual
dim(physeq2@tax_table) #6412 ASVs pre-rarefying
head(physeq2@otu_table)

batch1_rawasv<-colSums(physeq2@otu_table != 0) #before rarefying count non-zero rows in each column
summary(batch1_rawasv) #median 36, mean 53
sd(batch1_rawasv) #65

#Jaccard distance
jac3.11<-distance(ps3.11, method='jaccard', binary=T)

#chao1 and faith pd
alpha_r$chao<-estimate_richness(ps3.11, measures=c('Chao1'))[1] #add to alpha dataframe
kruskal.test(alpha_r$chao$Chao1, sample_data(ps3.11)$Species) #Kruskal-Wallis chi-squared = 19.324, df = 13, p-value = 0.1134
library(btools) #for faiths pd
testpd<-estimate_pd(ps3.11)
alpha_r$PD<-testpd$PD
kruskal.test(alpha_r$PD, sample_data(ps3.11)$Species) #Kruskal-Wallis chi-squared = 18.764, df = 13, p-value = 0.1306

#look at jacc pcoa--looks much like bray curtis
ord_j<-ordinate(ps3.11, "PCoA", "jaccard")
plot_ordination(ps3.11, ord_j, color='Species', shape = "State", title = "", type='samples') +
  scale_color_manual(values=pal) + geom_point(size=3) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(shape="Locality", colour="Species")

#PERMANOVA on jaccard
adonis2(jac3.11 ~ Species, data=md3.11) #r2=0.09722, p=0.001****
plot(betadisper(jac3.11, md3.11$Species))
permutest(betadisper(jac3.11, md3.11$Species)) #F=2.4086,p=0.009***
adonis2(jac3.11 ~ State, data=md3.11) #r2=0.01438, p=0.001****
plot(betadisper(jac3.11, md3.11$State))
permutest(betadisper(jac3.11, md3.11$State)) #F=5.6067,p=0.001***
adonis2(jac3.11 ~ Year, data=md3.11) #r2=0.01693, p=0.001****
plot(betadisper(jac3.11, md3.11$Year))
permutest(betadisper(jac3.11, md3.11$Year)) #F=5.889,p=0.004***

#using by=margin
adonis2(bray3.11 ~ Species+State+Year, data=md3.11, by='margin')
#species r2=0.08961, p=0.048
#state r2=0.01390, p=0.002
#year r2=0.01246, p=0.385
adonis2(bray3.11 ~ Species*State, data=md3.11, by='margin')

adonis2(jac3.11 ~ Species+State+Year, data=md3.11, by='margin')
#species r2=0.09336, p=0.001
#state r2=0.01103, p=0.001
#year r2=0.01531, p=0.001

adonis2(uni3.11 ~ Species+State+Year, data=md3.11, by='margin')
#species r2=0.10333, p=0.001
#state r2=0.00915, p=0.014
#year r2=0.01545, p=0.034

adonis2(wuni3.11 ~ Species+State+Year, data=md3.11, by='margin')
#species r2=0.08452, p=0.380
#state r2=0.01748, p=0.039
#year r2=0.00776, p=0.732



#for supplement
pdf(file='betadisper_3_2017-2019_v2.pdf', height=5, width=5)
myplotbetadisper(betadisper(jac3.11, md3.11$Species), ellipse = TRUE, hull = FALSE, 
                 fillrect=col.fill.rect, seg.lwd = 0.5,
                 alphaPoints=0.7, labPoints=NULL, sub='',
                 main= "Jaccard", coltextrect=pal, label.cex=0.6)
dev.off()

adonis2(jac3.11 ~ Species+State+Year, data=md3.11, by='margin') #working right here...can report margin results for batch 1, but not full dataset?

#export distance matricies for mantel tests with cytb
#batch 1 dataset
saveRDS(bray3.11, file='bray_b1.rds')
saveRDS(jac3.11, file='jac_b1.rds')
saveRDS(uni3.11, file='unifrac_b1.rds')
saveRDS(wuni3.11, file='wunifrac_b1.rds')
#grab batch 1 NY subset
psny<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$State=='NY')], ps3.11)
#remove zero sum taxa
psny<-prune_taxa(taxa_sums(psny) > 0, psny)
#recalculate distance matrices
bray_ny<-distance(psny, method='bray')
jac_ny<-distance(psny, method='jaccard', binary=T)
uni_ny<-distance(psny, method='unifrac')
wuni_ny<-distance(psny, method='wunifrac')
#save
saveRDS(bray_ny, file='bray_ny.rds')
saveRDS(jac_ny, file='jac_ny.rds')
saveRDS(uni_ny, file='unifrac_ny.rds')
saveRDS(wuni_ny, file='wunifrac_ny.rds')

#see how many of the 14 species carry each core bacteria
table(ps3.11@sam_data$Species)
head(core30) #core w/prevalence > 30%
amre<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='AMRE')],ps3.11)
amre<-prune_taxa(taxa_sums(amre) > 0, amre) #remove zero sum taxa
amre_c<-as.data.frame(amre@tax_table)
amre_c<-amre_c[which(row.names(amre_c) %in% rownames(core30)),] 
dim(amre_c) #all 9

baww<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='BAWW')],ps3.11)
baww<-prune_taxa(taxa_sums(baww) > 0, baww) #remove zero sum taxa
baww_c<-as.data.frame(baww@tax_table)
baww_c<-baww_c[which(row.names(baww_c) %in% rownames(core30)),] #all 9
dim(baww_c) #all 9

blbw<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='BLBW')],ps3.11)
blbw<-prune_taxa(taxa_sums(blbw) > 0, blbw) #remove zero sum taxa
blbw_c<-as.data.frame(blbw@tax_table)
blbw_c<-blbw_c[which(row.names(blbw_c) %in% rownames(core30)),] #all 9
dim(blbw_c) #all 9

btbw<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='BTBW')],ps3.11)
btbw<-prune_taxa(taxa_sums(btbw) > 0, btbw) #remove zero sum taxa
btbw_c<-as.data.frame(btbw@tax_table)
btbw_c<-btbw_c[which(row.names(btbw_c) %in% rownames(core30)),] #all 9
dim(btbw_c) #all 9

btnw<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='BTNW')],ps3.11)
btnw<-prune_taxa(taxa_sums(btnw) > 0, btnw) #remove zero sum taxa
btnw_c<-as.data.frame(btnw@tax_table)
btnw_c<-btnw_c[which(row.names(btnw_c) %in% rownames(core30)),] #all 9
dim(btnw_c) #all 9

cawa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='CAWA')],ps3.11)
cawa<-prune_taxa(taxa_sums(cawa) > 0, cawa) #remove zero sum taxa
cawa_c<-as.data.frame(cawa@tax_table)
cawa_c<-cawa_c[which(row.names(cawa_c) %in% rownames(core30)),] #all 9
dim(cawa_c) #all 9

coye<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='COYE')],ps3.11)
coye<-prune_taxa(taxa_sums(coye) > 0, coye) #remove zero sum taxa
coye_c<-as.data.frame(coye@tax_table)
coye_c<-coye_c[which(row.names(coye_c) %in% rownames(core30)),] #all 9
dim(coye_c) #all 9

cswa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='CSWA')],ps3.11)
cswa<-prune_taxa(taxa_sums(cswa) > 0, cswa) #remove zero sum taxa
cswa_c<-as.data.frame(cswa@tax_table)
cswa_c<-cswa_c[which(row.names(cswa_c) %in% rownames(core30)),] #all 9
dim(cswa_c) #all 9

howa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='HOWA')],ps3.11)
howa<-prune_taxa(taxa_sums(howa) > 0, howa) #remove zero sum taxa
howa_c<-as.data.frame(howa@tax_table)
howa_c<-howa_c[which(row.names(howa_c) %in% rownames(core30)),] 
dim(howa_c) #only 7 of 9
rownames(core30)[(which(!(rownames(core30) %in% rownames(howa_c))))] #"5648dccee530d68ceb3e4d7d22cf8756"-pseudomonas, "90f44084459949dd4fd965e1fc355026"-"g__1174-901-12" 

mawa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='MAWA')],ps3.11)
mawa<-prune_taxa(taxa_sums(mawa) > 0, mawa) #remove zero sum taxa
mawa_c<-as.data.frame(mawa@tax_table)
mawa_c<-mawa_c[which(row.names(mawa_c) %in% rownames(core30)),] 
dim(mawa_c) #all 9

mywa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='MYWA')],ps3.11)
mywa<-prune_taxa(taxa_sums(mywa) > 0, mywa) #remove zero sum taxa
mywa_c<-as.data.frame(mywa@tax_table)
mywa_c<-mywa_c[which(row.names(mywa_c) %in% rownames(core30)),] 
dim(mywa_c) #all 9

nawa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='NAWA')],ps3.11)
nawa<-prune_taxa(taxa_sums(nawa) > 0, nawa) #remove zero sum taxa
nawa_c<-as.data.frame(nawa@tax_table)
nawa_c<-nawa_c[which(row.names(nawa_c) %in% rownames(core30)),] 
dim(nawa_c) #all 9

nopa<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='NOPA')],ps3.11)
nopa<-prune_taxa(taxa_sums(nopa) > 0, nopa) #remove zero sum taxa
nopa_c<-as.data.frame(nopa@tax_table)
nopa_c<-nopa_c[which(row.names(nopa_c) %in% rownames(core30)),] 
dim(nopa_c) #all 9

oven<-prune_samples(rownames(ps3.11@sam_data)[which(ps3.11@sam_data$Species=='OVEN')],ps3.11)
oven<-prune_taxa(taxa_sums(oven) > 0, oven) #remove zero sum taxa
oven_c<-as.data.frame(oven@tax_table)
oven_c<-oven_c[which(row.names(oven_c) %in% rownames(core30)),] 
dim(oven_c) #only 7 of 9
rownames(core30)[(which(!(rownames(core30) %in% rownames(oven_c))))]
#"9477d45e83a2f8d1f61fd3cb9db73a0f"-f__Beijerinckiaceae "90f44084459949dd4fd965e1fc355026"-"g__1174-901-12" 

#make individual-level relative abundance plots for supp?
level_order2<-data.frame(id=mglom3.11$Sample,species=mglom3.11$sample_Species)
orders<-data.frame(species=level_order, order=1:15)
orders<-orders[-14,] #get rid of wewa
merge <- left_join(orders, level_order2, by = "species")
merge$Sample<-merge$id
orders[14,2]<-14
mglom3.11$order<-match(mglom3.11$sample_Species, orders$species)
mglom3.11$order2<-paste(mglom3.11$order, mglom3.11$Sample, sep='_')
mglom3.11$new<-mglom3.11$order
mglom3.11$new <- car::recode(mglom3.11$new, '10=91;11=92;12=93;13=94;14=95')
table(mglom3.11$new)      
mglom3.11$order3<-paste(mglom3.11$new, mglom3.11$Sample, sep='_')

#but use merged version, duh not this one.
pdf(file='~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/figs/ind_abund_micro.pdf', height=5, width=12)
ggplot(mglom3.11, aes(x = order3, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('azure3','turquoise4',
                             'rosybrown3','thistle1','violetred4',
                             'rosybrown2','wheat4','sandybrown','black')) +
  ggtitle('16s Batch 1') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  theme(legend.position="none") +
  geom_vline(xintercept = 14.5) + geom_vline(xintercept=29.5) + geom_vline(xintercept=42.5) +
  geom_vline(xintercept =51.5) +geom_vline(xintercept =63.5) +geom_vline(xintercept =74.5) +
  geom_vline(xintercept =87.5) + geom_vline(xintercept =103.5) +geom_vline(xintercept =109.5) +
  geom_vline(xintercept =121.5) + geom_vline(xintercept =130.5) + geom_vline(xintercept =141.5) +
  geom_vline(xintercept =150.5) + geom_vline(xintercept =161.5)
dev.off()

table(md3.11$Species)
level_order

save.image('phyloseq4_revision.RData')
load('phyloseq4_revision.RData')
