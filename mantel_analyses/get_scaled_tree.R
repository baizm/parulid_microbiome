library(ape)
library(brranching)
library(phangorn)

#set wd
setwd("~/Documents/Toews_Lab/manuscripts/parulidGM/reviews:revision/eliot_example/")

tree_full <- read.tree("~/Documents/Toews_Lab/16s_merged/TreeCmp/newick_warblers.newick")

#replace codes with scientific names
addedNames <- read.csv("addedNames_full.csv")
addedNames$species <- sub(" ", "_", addedNames$species)
merged <- merge(addedNames, data.frame(order=1:length(tree_full$tip.label), code=tree_full$tip.label))
merged <- merged[order(merged$order),]
tree_full$tip.label <- merged$species

#add some node labels. label the divergence between ovenbird and everything else (10.1 mya from time tree)
tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(ovenbird.everything=c("Seiurus","Setophaga")))

#and between cardellina and setophaga (5.9 mya from time tree)
tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(cardellina.setophaga=c("Cardellina","Setophaga")))

#and between redstart and rest of setophaga (5.3)
#tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(redstart.setophaga=c("Setophaga_ruticilla","Setophaga_virens")))

#and between geothlypis and leiothlypis (7.54)
tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(geo.leio=c("Geothlypis","Leiothlypis")))

#and between virens and coronata (5.3)
tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(virens.coro=c("Setophaga_virens","Setophaga_coronata")))

#and between pensylvanica and fusca (3.87)
tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(pen.fus=c("Setophaga_pensylvanica","Setophaga_fusca")))

#and between magnolia and americana (5.6)
tree_full <- makeNodeLabel(tree_full, method="user", nodeList=list(mag.am=c("Setophaga_magnolia","Setophaga_americana")))

#now prep a data frame for bladj
agesDF <- data.frame(a=c("ovenbird.everything","cardellina.setophaga","redstart.setophaga","geo.leio",
                         "virens.coro","pen.fus","mag.am"),
                     b=c(10.1, 5.9, 5.3, 7.54, 5.3, 3.87, 5.6))

#run bladj
scaled_full <- rbladj(tree=tree_full, ages=agesDF)

#this drops the upper case from the names and doesn't make a perfectly ultrametric tree. fix those
scaled_full$tip.label <- tree_full$tip.label
scaled_full <- nnls.tree(cophenetic(scaled_full), scaled_full, rooted=TRUE)

write.tree(scaled_full, "scaledWarbler_full.tre")
