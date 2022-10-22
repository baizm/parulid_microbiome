setwd("~/Documents/Toews_Lab/16s_data/Jan2020/silva/phyloseq/decontaminated/qiime_diversity/core-metrics-results_filtered/grouped_beta")
library(ctc)

#read in grouped distance matrices
bray<-read.csv('exported_braycurtis/distance-matrix.tsv', sep='\t', header=T, row.names=1)
uni<-read.csv('exported_unifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
wuni<-read.csv('exported_weightedunifrac/distance-matrix.tsv', sep='\t', header=T, row.names=1)
jac<-read.csv('exported_jaccard/distance-matrix.tsv', sep='\t', header=T, row.names=1)

#transform to matrix
bray<-as.matrix(bray)
bray<-as.dist(bray, diag=T, upper=F) #dist object for hclust

uni<-as.matrix(uni)
uni<-as.dist(uni, diag=T, upper=F)

wuni<-as.matrix(wuni)
wuni<-as.dist(wuni, diag=T, upper=F)

jac<-as.matrix(jac)
jac<-as.dist(jac, diag=T, upper=F)

#create UPGMA dendrograms by clustering
bd<-hclust(bray, method='average')
ud<-hclust(uni, method='average')
wd<-hclust(wuni, method='average')
jd<-hclust(jac, method='average')
plot(bd, hang=-1)
plot(ud, hang=-1)
plot(wd, hang=-1)
plot(jd, hang=-1)

#write to newick format for TreeCmp
bn<-hc2Newick(bd, flat=T)
un<-hc2Newick(ud, flat=T)
wn<-hc2Newick(wd, flat=T)
jn<-hc2Newick(jd, flat=T)

write.table(bn, 'newick_bray_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(un, 'newick_unifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(wn, 'newick_weightedunifrac_grouped.newick', quote=F, col.names=F, row.names=F)
write.table(jn, 'newick_jaccard_grouped.newick', quote=F, col.names=F, row.names=F)
