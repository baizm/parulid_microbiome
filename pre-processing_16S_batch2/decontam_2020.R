setwd("~/Documents/Toews_Lab/16s_2020samples/qiime_workflow/decontam/")
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
otu_table<-read.csv('feature-table.tsv', sep='\t', header=T, skip=1, check.names=FALSE)
#rearrange columns so negatives are at end to match sampledata (don't know if necessary)
otu_table<-otu_table[,c(1,6:158,2:5)]
#convert to matrix
otumat<-as.matrix(otu_table[,2:158])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
###-------read in taxonomy table-----------------
tax_table<-read.csv('silva_taxonomy.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('../../metadata/metadata_16s_2020samples.tsv', sep='\t', header=T,)
rownames(sampledata)<-sampledata$id
sampledata<-sampledata[,2:25]
#change Type col
sampledata$Type<-c(rep('pos', 153), rep('neg', 4))
###-------read in nwk tree with ape package-----------
tree<-read.tree('tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
physeq = phyloseq(OTU, TAX, SAM, tree)
#sums of OTUs across samples
summary(taxa_sums(physeq))



####------remove decontaminant ASVs--------##########
#add an is.neg column to sampledata
sample_data(physeq)$is.neg <- grepl('neg', sample_data(physeq)$Type) 
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant) #203 contaminant ASVs
contams<-contamdf.prev[which(contamdf.prev$contaminant),]

#subset otu table to only contain 203 contaminant ASVs
otu_table_contams<-otu_table(physeq)[which(rownames(otu_table(physeq)) %in% rownames(contams)),]
#look at last 4 cols where negs are
otu_table_contams[,154:157] #no overlap between two PCR negs and PCR negs and extr negs

#Let’s take a look at the number of times several of these taxa were observed in negative controls and positive samples
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(grepl('neg', sample_data(physeq)$Type), ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Type == "pos", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#there is at least one "non-contaminant" that hos a high prev in neg but low prev in sample...find it?
df.pa[which(df.pa$pa.neg==1 & df.pa$pa.pos<5 & !df.pa$contaminant),] #they are all pos=0
tail(df.pa[order(df.pa$pa.neg, df.pa$pa.pos),])
df.pa[which(df.pa$contaminant==F & df.pa$pa.neg==1),] #many are in neg, not in pos, but still not ID as contamination
#take a look
otu_table(physeq)[which(rownames(otu_table(physeq))=='3e237d87fa880fc3a21be7a0445e1bec'),]
#so its all ps.neg=1 and pa.pos=0 and contaminnant=F, but should be true. p=NA bc they aren't observed in the sample

##compile list on contaminant ASVs
#here are ASV's that should also be removed bc they were only observed in negatives
neg_only<-rownames(df.pa[which(df.pa$pa.neg>0 & df.pa$pa.pos==0 & !df.pa$contaminant),]) #156 ASVs
#here are ASVs identified by isContaminant as contaminants (threshold=0.5)
contam_ids<-rownames(contams) #203 ASVs
#compile vector of ASV ids to remove
asv_to_remove<-c(neg_only, contam_ids) #all contaminant ASVs, 359 of them
#for qiime2 to filter
asv_to_remove<-data.frame(featureid=c(asv_to_remove)) #has to have this colname for qiime2
write.csv(asv_to_remove,'asv_to_remove.txt', quote=F, sep='\t', row.names=F, col.names=T)
