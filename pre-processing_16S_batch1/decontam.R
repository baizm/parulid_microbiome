setwd("~/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq")
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

########-------create physeq object---------------------########
####------read in OTU table, has to be a matrix---------------------
otu_table<-read.csv('../feature-table.tsv', sep='\t', header=T, skip=1)
#rearrange columns so negatives are at end to match sampledata (don't know if necessary)
otu_table<-otu_table[,c(1,6:256,2:5)]
#convert to matrix
otumat<-as.matrix(otu_table[,2:256])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
#correct negative names so they match other tables
colnames(otumat)[252:255]<-c('1-neg-PCR', '2-neg-PCR', '7-neg-extr', '8-neg-extr')
###-------read in taxonomy table-----------------
tax_table<-read.csv('../silva_taxonomy.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#it wasn't filtered, so i'll have to delete mt,cp,un...
tax_table2<-tax_table2[which(tax_table2$Feature.ID %in% otu_table$X.OTU.ID),]
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('../../metadata2_16s_ny.tsv', sep='\t', header=T,)
rownames(sampledata)<-sampledata$id
sampledata<-sampledata[,2:29]
###-------read in nwk tree with ape package-----------
tree<-read.tree('../tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
physeq = phyloseq(OTU, TAX, SAM, tree)
#sums of OTUs across samples
taxa_sums(physeq)
#alpha - tells me I don't have singletons...is this suspicious?
plot_richness(physeq, x="Species", color="Species",measures="Shannon") + theme(legend.position="none") 
sums<-sample_sums(physeq)
physeq@sam_data$depth<-sums

#for molecular ecology course 2022
sums2<-data.frame(row.names=rownames(physeq@sam_data), species=physeq@sam_data$Species, depth=physeq@sam_data$depth, site=physeq@sam_data$State)
fd<-c('AMRE', 'BTNW', 'CAWA')
dave<-physeq@sam_data[physeq@sam_data$Year=='2019']
dave<-dave[dave$Species %in% fd,]
dave<-dave[dave$depth > 3000,]

write.csv(rownames(dave), '~/Documents/Toews_Lab/MolecEcol_course/dave_samples.csv', quote=F)
write.csv(dave[,1:28], '~/Documents/Toews_Lab/MolecEcol_course/dave_metadata.csv', quote=F)

fd<-c('AMRE', 'BTNW')
dave<-physeq@sam_data[physeq@sam_data$Year=='2019']
dave<-dave[dave$Species %in% fd,]
dave[order(dave$depth, decreasing=T)]
write.csv(rownames(dave), '~/Documents/Toews_Lab/MolecEcol_course/dave_samples3.csv', quote=F)
write.csv(dave[,1:28], '~/Documents/Toews_Lab/MolecEcol_course/dave_metadata3.csv', quote=F)
####------remove decontaminant ASVs--------##########
#add an is.neg column to sampledata
sample_data(physeq)$is.neg <- grepl('neg', sample_data(physeq)$Type) 
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant) #52 contaminant ASVs
contams<-contamdf.prev[which(contamdf.prev$contaminant),]

#subset otu table to only contain 52 contaminant ASVs
otu_table_contams<-otu_table(physeq)[which(rownames(otu_table(physeq)) %in% rownames(contams)),]
#look at last 4 rows where negs are
otu_table_contams[,252:255] #no overlap between two PCR negs and PCR negs and extr negs

#Let’s take a look at the number of times several of these taxa were observed in negative controls and positive samples
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(grepl('neg', sample_data(physeq)$Type), ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Type == "sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#there is at least one "non-contaminant" that hos a high prev in neg but low prev in sample...find it?
which(df.pa$pa.neg==1 & df.pa$pa.pos<5 & !df.pa$contaminant) #they are all pos=0
tail(df.pa[order(df.pa$pa.neg, df.pa$pa.pos),])
df.pa[which(df.pa$contaminant==F & df.pa$pa.neg==1),] #many are in neg, not in pos, but still not ID as contamination
#take a look
otu_table(physeq)[which(rownames(otu_table(physeq))=='2e2254ef935ea621fc582ebc68576bb7'),]
#so its all ps.neg=1 and pa.pos=0 and contaminnant=F, but should be true. p=NA bc they aren't observed in the sample i think?

##compile list on contaminant ASVs
#here are ASV's that should also be removed bc they were only observed in negatives
neg_only<-rownames(df.pa[which(df.pa$pa.neg>0 & df.pa$pa.pos==0 & !df.pa$contaminant),])
#here are ASVs identified by isContaminant as contaminants (threshold=0.5)
contam_ids<-rownames(contams)
#compile vector of ASV ids to remove
asv_to_remove<-c(neg_only, contam_ids) #all contaminant ASVs, 87 of them
#for qiime2 to filter
asv_to_remove<-data.frame(feature_id=c(asv_to_remove)) #has to have this colname for qiime2
write.csv(asv_to_remove,'asv_to_remove.txt', quote=F, sep='\t', row.names=F, col.names=T)

 #remove 87 contaminant ASVs from physeq object
#make logical vector to prune from ("a logical vector where the kept taxa are TRUE, and length is equal to the number of taxa in object x")
v<-which(!(taxa_names(physeq) %in% asv_to_remove))
v2<-taxa_names(physeq)[c(v)]
#make sure none of the 89 ASVs are in v2
which(asv_to_remove %in% v2) #none, it worked!
ps<-prune_taxa(v2, physeq)

