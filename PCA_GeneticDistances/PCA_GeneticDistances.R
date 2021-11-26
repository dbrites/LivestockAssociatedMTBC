---
title: "livestock MTBC nomenclature- PCA and genetic distances"
author: "Daniela Brites"
date: "4/9/2021"
---
  
#R/3.6.0
library(adegenet)
library(ape)
library(sqldf)
library(glue)
library(dplyr)

##PCA
alignmentMTBC_genlight <- fasta2genlight("genomes_mtbc_1226_WGS.fasta_var_N.fasta",chunk=50, parallel=FALSE,saveNbAlleles=T)
#remove outgroup and remove sequences identical as identified by raxml
length(alignmentMTBC_genlight$ind.names)
#[1] 1227

alignmentMTBC_genlight1<- alignmentMTBC_genlight[indNames(alignmentMTBC_genlight) != "G00157"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G00205"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G37371"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G00240"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G00220"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G47573"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G26462"]
alignmentMTBC_genlight1<- alignmentMTBC_genlight1[indNames(alignmentMTBC_genlight1) != "G60801"]
length(alignmentMTBC_genlight1$ind.names)

#PCA of the complex
pca_mtbc <- glPca(alignmentMTBC_genlight1,parallel=FALSE,nf=20) 

pca_mtbc.scores <- as.data.frame(pca_mtbc$scores)
rownames(pca_mtbc.scores)
names_pca <-as.data.frame(rownames(pca_mtbc.scores),stringsAsFactors =F)
dim(names_pca)
colnames(names_pca) <- "g_number"

#read metadata
metadata <-read.csv(file="metadata1226Fig1.txt",sep="\t",header=T)
#order metadata gnumbers as in pca
metadata_ordered<- joined.data[order(match(joined.data$g_number,names_pca$g_number)), ]

#attribute a population to each of the gnumbers in the PCA
pca_mtbc.scores1and2 <-pca_mtbc.scores[,c(1,2)]
pca_mtbc.scores1and2$GROUP <-metadata_ordered$lineage
dim(metadata_ordered)

pca_mtbc$eig[1]/sum(pca_mtbc$eig)*100
pca_mtbc$eig[2]/sum(pca_mtbc$eig)*100
pca_mtbc$eig[3]/sum(pca_mtbc$eig)*100
barplot(pca_mtbc$eig/sum(pca_mtbc$eig)*100, main="eigenvalues", col=heat.colors(length(pca_mtbc$eig)))

#Plot PCA

d <- c( "Mbovis"="#EFA0CF",  "L1" = "#ff00ff", "L2" = "#0000ff", "L3" = "#a000cc",  "L4" = "#ff0000", "L5" = "#663200", "L6" = "#00cc33", "L7" = "#ede72e",  "6AChimp" = "gray77", "D" = "#000000","L9" =  "#006400", "L8" ="#FF7F00", "C"= "grey","Morygis"="lightpink2", "D"="grey","pinnipedii"="grey","D"="#000000","Mcaprae" = "cyan","Mbovis PZA sus"="purple")

library(ggplot2)
set.seed(9)

d_grey <- c( "Mbovis"="grey27",  "L1" = "#ff00ff", "L2" = "#0000ff", "L3" = "#a000cc",  "L4" = "#ff0000", "L5" = "#663200", "L6" = "#00cc33", "L7" = "#ede72e",  "6AChimp" = "gray77", "D" = "gray77","L9" =  "#006400", "L8" ="#FF7F00", "C"= "gray77","Morygis"="black", "D"="gray77","pinnipedii"="grey","D"="gray77","Mcaprae" = "gray87","Mbovis PZA sus"="gray50")
p <- p+geom_point(size=3,alpha=0.5)+geom_hline(yintercept=0,size=0.02,color="grey")+geom_vline(xintercept=0,size=0.02,color="grey")+xlab("PC1 (24.4%)") + ylab("PC2 (11.3%)")
p2 <- p+geom_point(size=3,alpha=0.5)+geom_hline(yintercept=0,size=0.02,color="grey")+geom_vline(xintercept=0,size=0.02,color="grey")+xlab("PC1 (24.4%)") + ylab("PC2 (11.3%)")
p3 <-p2+ scale_color_manual(values = d_grey) + theme(legend.position = "none")



#Calculate and plot pair-wise genetic distances using the 1226 dataset
distdna_1126MTBC_pairwisedel<- dist.dna(alignmentMTBC_genlight1,model="N",variance=F,as.matrix=T,pairwise.deletion =T)

capraeGnumbers <-as.character(metadata_ordered[metadata_ordered$lineage=="Mcaprae",1])
pzaSusGnumbers <-as.character(metadata_ordered[metadata_ordered$lineage=="Mbovis PZA sus",1])
bovisGnumbers <- as.character(metadata_ordered[metadata_ordered$lineage=="Mbovis",1])
orygisGnumbers <- as.character(metadata_ordered[metadata_ordered$lineage=="Morygis",1])

dist_capraeVSpzaSus<-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% capraeGnumbers, rownames(distdna_1126MTBC_pairwisedel) %in% pzaSusGnumbers]
mean(dist_capraeVSpzaSus)

dist_capraeVbovis<-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% capraeGnumbers , rownames(distdna_1126MTBC_pairwisedel) %in% bovisGnumbers]
mean(dist_capraeVbovis)

dist_bovisVSorygis <-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% orygisGnumbers, rownames(distdna_1126MTBC_pairwisedel) %in% bovisGnumbers]
mean(dist_bovisVSorygis)

dist_bovisVSpzaSus <-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% pzaSusGnumbers , rownames(distdna_1126MTBC_pairwisedel) %in% bovisGnumbers]
mean(dist_bovisVSpzaSus)

#build data frame in order to plot the distances
matrixDist=c(dist_bovisVSorygis,dist_capraeVbovis,dist_bovisVSpzaSus,dist_capraeVSpzaSus)
matrixName=c(rep("Mbovis-Morygis",length(dist_bovisVSorygis)),rep("Mbovis-Mcaprae",length(dist_capraeVbovis)),rep("Mbovis(PZAres)-Mbovis(PZAsus)",length(dist_bovisVSpzaSus)),rep("Mcaprae-Mbovis(PZAsus)",length(dist_capraeVSpzaSus)))
table1 <-data.frame(distances=length(matrixDist),Species=1,stringsAsFactors = FALSE)
for (dist in 1:length(matrixDist)){
  table1[dist,1] =matrixDist[dist]
  table1[dist,2] =matrixName[dist]
}

#do the same for L1, L2, L3, L4,L5,L6
L2Gnumbers <-as.character(joined.data[joined.data$lineage=="L2",1])
L3Gnumbers <-as.character(joined.data[joined.data$lineage=="L3",1])
L4Gnumbers <- as.character(joined.data[joined.data$lineage=="L4",1])
L1Gnumbers <- as.character(joined.data[joined.data$lineage=="L1",1])
L5Gnumbers <- as.character(joined.data[joined.data$lineage=="L5",1])
L6Gnumbers <- as.character(joined.data[joined.data$lineage=="L6",1])

dist_L2L3<-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% L2Gnumbers, rownames(distdna_1126MTBC_pairwisedel) %in% L3Gnumbers]
mean(dist_L2L3)

dist_L2L4<-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% L2Gnumbers , rownames(distdna_1126MTBC_pairwisedel) %in% L4Gnumbers]
mean(dist_L2L4)

dist_L4L1 <-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% L1Gnumbers, rownames(distdna_1126MTBC_pairwisedel) %in% L4Gnumbers]
mean(dist_L4L1)

dist_L4L3 <-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% L3Gnumbers , rownames(distdna_1126MTBC_pairwisedel) %in% L4Gnumbers]
mean(dist_L4L3)

dist_L5L6 <-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% L5Gnumbers , rownames(distdna_1126MTBC_pairwisedel) %in% L6Gnumbers]
mean(dist_L5L6)

matrixDist2=c(dist_L4L1,dist_L2L4,dist_L4L3,dist_L2L3,dist_L5L6)
matrixName2=c(rep("L4-L1",length(dist_L4L1)),rep("L4-L2",length(dist_L2L4)),rep("L4-L3",length(dist_L4L3)),rep("L2-L3",length(dist_L2L3)),rep("L5-L6",length(dist_L5L6)))
table2 <-data.frame(distances=length(matrixDist2),Species=1,stringsAsFactors = FALSE)
for (dist in 1:length(matrixDist2)){
  table2[dist,1] =matrixDist2[dist]
  table2[dist,2] =matrixName2[dist]
}

#combine both tables
table3<-(rbind(table1,table2))

gnumbersProtoBeijing <-c("G62034","G00056","G00156","G00828","G00910","G00912","G00917","G02404","G03153","G03836","G06122","G06123","G06124","G06125","G06126","G06127","G06128","G06129","G06130","G06131","G06132","G06133","G06134","G06135","G06136","G06137")
dist_L2Protobeijing<-distdna_1126MTBC_pairwisedel[rownames(distdna_1126MTBC_pairwisedel) %in% L2Gnumbers, rownames(distdna_1126MTBC_pairwisedel) %in% gnumbersProtoBeijing]
hist(dist_L2Protobeijing)

which(rownames(dist_L2Protobeijing) %in% gnumbersProtoBeijing)
#[1]  62 105 164

#remove those from the rows
distL2noProto.Protobeijing <-dist_L2Protobeijing[-c(62,105,164),]
dim(dist_L2Protobeijing)
dim(distL2noProto.Protobeijing)
hist(distL2noProto.Protobeijing)

matrixDist3=c(distL2noProto.Protobeijing)
matrixName3=c(rep("L2-ProtoBeijing",length(distL2noProto.Protobeijing)))
table4 <-data.frame(distances=length(matrixDist3),Species=1,stringsAsFactors = FALSE)

for (dist in 1:length(matrixDist3)){
  table4[dist,1] =matrixDist3[dist]
  table4[dist,2] =matrixName3[dist]
}

#combine both tables
table5<-(rbind(table1,table2,table4))

#Plot distances
plot_dist <-ggplot(table5, aes(x = distances, y=Species)) +
  geom_density_ridges(aes(fill=Species),)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "SNP pair-wise distances")+ 
  theme(axis.text.x = element_text(face="bold",size=12),axis.text.y = element_text(face="bold",size=12))



##Calculate genetic distances within La1
alignment749<- read.dna("genomes_boviscaprae_748_WGSvar_editedX.fasta",format="fasta")
#remove outgroup
alignment_748<-alignment749[rownames(alignment749) != "G08466",]
dist748_pairwisedel<- dist.dna(alignment_748,model="N",variance=F,as.matrix=T,pairwise.deletion =T)
hist(dist748_pairwisedel)
dist748_pairwisedel

metadata <-read.csv(file="TableS1.V1.txt",header=T,sep="\t",check.names=FALSE,stringsAsFactors = F)
dim(metadata)
#[1] 840  10

#remove L6
which(metadata$Clonal_Complex=="L6")
#[1] 466
metadata <- metadata[-466,]
dim(metadata)
#[1] 839  10

#remove orygis from metadata
table(metadata$Clonal_Complex)
which(metadata$Clonal_Complex=="orygis")
metadata748 <-metadata[-which(metadata$Clonal_Complex=="orygis"),]
dim(metadata748)




af2Gnumbers <-as.character(metadata748[metadata748$Clonal_Complex=="af2",5])
un2Gnumbers <-as.character(metadata748[metadata748$Clonal_Complex=="unk2",5])
un3Gnumbers <-as.character(metadata748[metadata748$Clonal_Complex=="unk3",5])
af1Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="af1",5])
un9Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="unk9",5])
eu2Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="eu2",5])
un4Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="unk4",5])
un5Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="unk5",5])
un7Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="unk7",5])
un6Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="unk6",5])
eu1Gnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="eu1",5])
pzaGnumbers <- as.character(metadata748[metadata748$Clonal_Complex=="PZA_sus",5])


distPzaAf2 <- dist748_pairwisedel[rownames(dist748_pairwisedel) %in% pzaGnumbers, rownames(dist748_pairwisedel) %in% af2Gnumbers]
distPza <- dist748_pairwisedel[rownames(dist748_pairwisedel) %in% pzaGnumbers, rownames(dist748_pairwisedel) %in% pzaGnumbers]
distAf2un2<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% af2Gnumbers, rownames(dist748_pairwisedel) %in% un2Gnumbers]
distAf2<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% af2Gnumbers, rownames(dist748_pairwisedel) %in% af2Gnumbers]
distUn2<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un2Gnumbers, rownames(dist748_pairwisedel) %in% un2Gnumbers]
distAf1un9<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% af1Gnumbers, rownames(dist748_pairwisedel) %in% un9Gnumbers]
distUn3un9<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un3Gnumbers, rownames(dist748_pairwisedel) %in% un9Gnumbers]
distUn3Af1<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un3Gnumbers, rownames(dist748_pairwisedel) %in% af1Gnumbers]
distUn3 <-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un3Gnumbers, rownames(dist748_pairwisedel) %in% un3Gnumbers]
distAf1<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% af1Gnumbers, rownames(dist748_pairwisedel) %in% af1Gnumbers]
distUn9<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un9Gnumbers, rownames(dist748_pairwisedel) %in% un9Gnumbers]
distEu2un4<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% eu2Gnumbers, rownames(dist748_pairwisedel) %in% un4Gnumbers]
distUn5Un4<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un5Gnumbers, rownames(dist748_pairwisedel) %in% un4Gnumbers]
distUn5eu2 <-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un5Gnumbers, rownames(dist748_pairwisedel) %in% eu2Gnumbers]
distEu2<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% eu2Gnumbers, rownames(dist748_pairwisedel) %in% eu2Gnumbers]
distUn4<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un4Gnumbers, rownames(dist748_pairwisedel) %in% un4Gnumbers]
distUn5<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un5Gnumbers, rownames(dist748_pairwisedel) %in% un5Gnumbers]
distUn7eu1 <-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un7Gnumbers, rownames(dist748_pairwisedel) %in% eu1Gnumbers]
distUn7un6 <-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un7Gnumbers, rownames(dist748_pairwisedel) %in% un6Gnumbers]
distUn6eu1 <-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un6Gnumbers, rownames(dist748_pairwisedel) %in% eu1Gnumbers]
distEu1<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% eu1Gnumbers, rownames(dist748_pairwisedel) %in% eu1Gnumbers]
distUn6<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un6Gnumbers, rownames(dist748_pairwisedel) %in% un6Gnumbers]
distUn7<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un7Gnumbers, rownames(dist748_pairwisedel) %in% un7Gnumbers]

matrixDistAll=c(distPzaAf2,distPza,distAf2un2,distAf2,distUn2,distAf1un9,distUn3un9,distUn3Af1,distUn3,distAf1,distUn9,distEu2un4,distUn5Un4,distUn5eu2,distEu2,distUn4,distUn5,distUn7eu1,distUn7un6,distUn6eu1,distEu1,distUn6,distUn7)

matrixNameAll <-c(rep("PZAsus-Af2",length(distPzaAf2)),
                  rep("Pza",length(distPza)),
                  rep("Af2-Unk2",length(distAf2un2)),
                  rep("Af2",length(distAf2)),
                  rep("Unk2",length(distUn2)),
                  rep("Af1-Unk9",length(distAf1un9)),
                  rep("Unk3-Unk9",length(distUn3un9)),
                  rep("Unk3-Af1",length(distUn3Af1)),
                  rep("Unk3",length(distUn3)),
                  rep("Af1",length(distAf1)),
                  rep("Unk9",length(distUn9)),
                  rep("Eu2-Unk4",length(distEu2un4)),
                  rep("Unk5-Unk4",length(distUn5Un4)),
                  rep("Unk5-Eu2",length(distUn5eu2)),
                  rep("Eu2",length(distEu2)),
                  rep("Unk4",length(distUn4)),
                  rep("Unk5",length(distUn5)),
                  rep("Eu1-Unk7",length(distUn7eu1)),
                  rep("Unk7-Unk6",length(distUn7un6)),
                  rep("Eu1-Unk6",length(distUn6eu1)),
                  rep("Eu1",length(distEu1)),
                  rep("Unk6",length(distUn6)),
                  rep("Unk7",length(distUn7)))
table <-data.frame(distances=length(matrixDistAll),Species=1,stringsAsFactors = FALSE)
for (dist in 1:length(matrixDistAll)){
  table[dist,1] =matrixDistAll[dist]
  table[dist,2] =matrixNameAll[dist]
}

library(ggplot2)
library(ggridges)

table$Species <-factor(table$Species,levels=c("PZAsus-Af2",
                                              "Af2-Unk2",
                                              "Af1-Unk9",
                                              "Unk3-Unk9",
                                              "Unk3-Af1",
                                              "Eu2-Unk4",
                                              "Unk5-Unk4",
                                              "Unk5-Eu2",
                                              "Eu1-Unk7",
                                              "Unk7-Unk6",
                                              "Eu1-Unk6",
                                              "Pza",
                                              "Af2",
                                              "Unk2",
                                              "Unk3",
                                              "Af1",
                                              "Unk9",
                                              "Eu2",
                                              "Unk4",
                                              "Unk5",
                                              "Eu1",
                                              "Unk6",
                                              "Unk7"
))

ggplot(table, aes(x = distances, y = Species)) +
  geom_density_ridges(aes(fill = Species))

ggplot(table, aes(x = distances, y=Species)) +
  geom_density_ridges(aes(fill=Species),)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "SNP pair-wise distances")+
  theme(axis.text.x = element_text(face="bold",size=12),axis.text.y = element_text(face="bold",size=12))


##Create second factor to order the labels.
table$Species2 <-factor(table$Species,levels=c("PZAsus-Af2",
                                               "Pza",
                                               "Af2-Unk2",
                                               "Af2",
                                               "Unk2",
                                               "Unk3-Af1",
                                               "Unk3-Unk9",
                                               "Af1-Unk9",
                                               "Unk3",
                                               "Af1",
                                               "Unk9",
                                               "Eu2-Unk4",
                                               "Unk5-Eu2",
                                               "Unk5-Unk4",
                                               "Eu2",
                                               "Unk4",
                                               "Unk5",
                                               "Eu1-Unk7",
                                               "Eu1-Unk6",
                                               "Unk7-Unk6",
                                               "Eu1",
                                               "Unk6",
                                               "Unk7"
))
ggplot(table, aes(x = distances, y = Species2)) +
  geom_density_ridges(aes(fill = Species2))

#put un5 and un5 and Eu2 together in one group
#put un6 and un7 and Eu1 together in another group




un5_4_eu2Gnumbers <-c(un4Gnumbers,un4Gnumbers,eu2Gnumbers)
un6_7_eu1Gnumbers <-c(un6Gnumbers,un7Gnumbers,eu1Gnumbers)

distUn4un5Eu2_Un6un7Eu1<-dist748_pairwisedel[rownames(dist748_pairwisedel) %in% un5_4_eu2Gnumbers, rownames(dist748_pairwisedel) %in% un6_7_eu1Gnumbers]
matrixDistAll_2=c(distPzaAf2,distPza,distAf2un2,distAf2,distUn2,distAf1un9,distUn3un9,distUn3Af1,distUn3,distAf1,distUn9,distUn4un5Eu2_Un6un7Eu1,distEu2,distUn4,distUn5,distEu1,distUn6,distUn7)
matrixNameAll_2 <-c(rep("PZAsus-Af2",length(distPzaAf2)),
                    rep("Pza",length(distPza)),
                    rep("Af2-Unk2",length(distAf2un2)),
                    rep("Af2",length(distAf2)),
                    rep("Unk2",length(distUn2)),
                    rep("Af1-Unk9",length(distAf1un9)),
                    rep("Unk3-Unk9",length(distUn3un9)),
                    rep("Unk3-Af1",length(distUn3Af1)),
                    rep("Unk3",length(distUn3)),
                    rep("Af1",length(distAf1)),
                    rep("Unk9",length(distUn9)),
                    rep("Unk4unk5Eu2-Unk6unk7Eu1",length(distUn4un5Eu2_Un6un7Eu1)),
                    rep("Eu2",length(distEu2)),
                    rep("Unk4",length(distUn4)),
                    rep("Unk5",length(distUn5)),
                    rep("Eu1",length(distEu1)),
                    rep("Unk6",length(distUn6)),
                    rep("Unk7",length(distUn7)))
table_2 <-data.frame(distances=length(matrixDistAll_2),Species=1,stringsAsFactors = FALSE)
for (dist in 1:length(matrixDistAll_2)){
  table_2[dist,1] =matrixDistAll_2[dist]
  table_2[dist,2] =matrixNameAll_2[dist]
}

table_2$Species3 <-factor(table_2$Species,levels=c("PZAsus-Af2",
                                                   "Pza",
                                                   "Af2-Unk2",
                                                   "Af2",
                                                   "Unk2",
                                                   "Unk3-Af1",
                                                   "Unk3-Unk9",
                                                   "Af1-Unk9",
                                                   "Unk3",
                                                   "Af1",
                                                   "Unk9",
                                                   "Unk4unk5Eu2-Unk6unk7Eu1",
                                                   "Eu2",
                                                   "Unk4",
                                                   "Unk5",
                                                   "Eu1",
                                                   "Unk6",
                                                   "Unk7"
))
plot <-ggplot(table_2, aes(x = distances, y = Species3)) +
  geom_density_ridges(aes(fill = Species3))

#Following reviews a few duplicated sequences were identified withih La1. After removing duplicates genetic distances were re-calculated and plotted. 

library(ape) 
alignment738<- read.dna("boviscaprae_738_noduplicates_WGS.fasta_var_editedX.fasta",format="fasta")
#remove outgroup
alignment_737<-alignment738[rownames(alignment738) != "G08466",] dist737_pairwisedel<- dist.dna(alignment_737,model="N",variance=F,as.matrix=T,pairwise.deletion =T) 
hist(dist737_pairwisedel) 
metadata <-read.csv(file="Table1_corrected.txt",header=T,sep="\t",check.names=FALSE,stringsAsFactors = F) dim(metadata)
#[1] 830 13 remove L6
which(metadata$Sublineage=="L6")
#[1] 462
metadata <- metadata[-462,] dim(metadata)
#[1] 829 12 remove orygis from metadata
table(metadata$Sublineage) which(metadata$Sublineage=="La3") metadataLa1La2 <-metadata[-which(metadata$Sublineage=="La3"),] dim(metadataLa1La2)
#[1] 738 13
af2Gnumbers <-as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="Af2",5]) 
un2Gnumbers <-as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown2",5]) 
un3Gnumbers <-as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown3",5]) 
af1Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="Af1",5]) 
un9Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown9",5]) 
eu2Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="Eu2",5]) 
un4Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown4",5]) 
un5Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown5",5]) 
un7Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown7",5]) 
un6Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="unknown6",5]) 
eu1Gnumbers <- as.character(metadataLa1La2[metadataLa1La2$ClonalComplex_Loiseau2020=="Eu1",5]) 
pzaGnumbers <- as.character(metadataLa1La2[metadataLa1La2$Sublineage=="La1.1",5]) 
distPzaAf2 <- dist737_pairwisedel[rownames(dist737_pairwisedel) %in% pzaGnumbers, rownames(dist737_pairwisedel) %in% af2Gnumbers] distPza <- dist737_pairwisedel[rownames(dist737_pairwisedel) %in% pzaGnumbers, rownames(dist737_pairwisedel) %in% pzaGnumbers] 
distAf2un2<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% af2Gnumbers, rownames(dist737_pairwisedel) %in% un2Gnumbers] 
distAf2<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% af2Gnumbers, rownames(dist737_pairwisedel) %in% af2Gnumbers] 
distUn2<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un2Gnumbers, rownames(dist737_pairwisedel) %in% un2Gnumbers] 
distAf1un9<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% af1Gnumbers, rownames(dist737_pairwisedel) %in% un9Gnumbers] 
distUn3un9<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un3Gnumbers, rownames(dist737_pairwisedel) %in% un9Gnumbers] 
distUn3Af1<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un3Gnumbers, rownames(dist737_pairwisedel) %in% af1Gnumbers] 
distUn3 <-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un3Gnumbers, rownames(dist737_pairwisedel) %in% un3Gnumbers] 
distAf1<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% 
af1Gnumbers, rownames(dist737_pairwisedel) %in% af1Gnumbers] 
distUn9<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un9Gnumbers, rownames(dist737_pairwisedel) %in% un9Gnumbers] 
distEu2un4<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% eu2Gnumbers, rownames(dist737_pairwisedel) %in% un4Gnumbers] 
distUn5Un4<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un5Gnumbers, rownames(dist737_pairwisedel) %in% un4Gnumbers] 
distUn5eu2 <-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un5Gnumbers, rownames(dist737_pairwisedel) %in% eu2Gnumbers] 
distEu2<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% eu2Gnumbers, rownames(dist737_pairwisedel) %in% eu2Gnumbers] 
distUn4<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un4Gnumbers, rownames(dist737_pairwisedel) %in% un4Gnumbers] 
distUn5<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un5Gnumbers, rownames(dist737_pairwisedel) %in% un5Gnumbers] 
distUn7eu1 <-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un7Gnumbers, rownames(dist737_pairwisedel) %in% eu1Gnumbers] 
distUn7un6 <-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un7Gnumbers, rownames(dist737_pairwisedel) %in% un6Gnumbers] 
distUn6eu1 <-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un6Gnumbers, rownames(dist737_pairwisedel) %in% eu1Gnumbers] 
distEu1<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% eu1Gnumbers, rownames(dist737_pairwisedel) %in% eu1Gnumbers] 
distUn6<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un6Gnumbers, rownames(dist737_pairwisedel) %in% un6Gnumbers] 
distUn7<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un7Gnumbers, rownames(dist737_pairwisedel) %in% un7Gnumbers] 
matrixDistAll=c(distPzaAf2,distPza,distAf2un2,distAf2,distUn2,distAf1un9,distUn3un9,distUn3Af1,distUn3,distAf1,distUn9,distEu2un4,distUn5Un4,distUn5eu2,distEu2,distUn4,distUn5,distUn7eu1,distUn7un6,distUn6eu1,distEu1,distUn6,distUn7) 
matrixNameAll <-c(rep("PZAsus-Af2",length(distPzaAf2)),
                rep("Pza",length(distPza)),
                rep("Af2-Unk2",length(distAf2un2)),
                rep("Af2",length(distAf2)),
                rep("Unk2",length(distUn2)),
                rep("Af1-Unk9",length(distAf1un9)),
                rep("Unk3-Unk9",length(distUn3un9)),
                rep("Unk3-Af1",length(distUn3Af1)),
                rep("Unk3",length(distUn3)),
                rep("Af1",length(distAf1)),
                rep("Unk9",length(distUn9)),
                rep("Eu2-Unk4",length(distEu2un4)),
                rep("Unk5-Unk4",length(distUn5Un4)),
                rep("Unk5-Eu2",length(distUn5eu2)),
                rep("Eu2",length(distEu2)),
                rep("Unk4",length(distUn4)),
                rep("Unk5",length(distUn5)),
                rep("Eu1-Unk7",length(distUn7eu1)),
                rep("Unk7-Unk6",length(distUn7un6)),
                rep("Eu1-Unk6",length(distUn6eu1)),
                rep("Eu1",length(distEu1)),
                rep("Unk6",length(distUn6)),
                rep("Unk7",length(distUn7))) table <-data.frame(distances=length(matrixDistAll),Species=1,stringsAsFactors = FALSE) for (dist in 1:length(matrixDistAll)){
  table[dist,1] =matrixDistAll[dist]
  table[dist,2] =matrixNameAll[dist]
}
library(ggplot2) library(ggridges) table$Species <-factor(table$Species,levels=c("PZAsus-Af2",
                                              "Af2-Unk2",
                                              "Af1-Unk9",
                                              "Unk3-Unk9",
                                              "Unk3-Af1",
                                              "Eu2-Unk4",
                                              "Unk5-Unk4",
                                              "Unk5-Eu2",
                                              "Eu1-Unk7",
                                              "Unk7-Unk6",
                                              "Eu1-Unk6",
                                              "Pza",
                                              "Af2",
                                              "Unk2",
                                              "Unk3",
                                              "Af1",
                                              "Unk9",
                                              "Eu2",
                                              "Unk4",
                                              "Unk5",
                                              "Eu1",
                                              "Unk6",
                                              "Unk7" )) ggplot(table, aes(x = distances, y = Species)) +
  geom_density_ridges(aes(fill = Species)) ggplot(table, aes(x = distances, y=Species)) +
  geom_density_ridges(aes(fill=Species),)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "SNP pair-wise distances")+
  theme(axis.text.x = element_text(face="bold",size=12),axis.text.y = element_text(face="bold",size=12))
##Create second factor to order the labels.
table$Species2 <-factor(table$Species,levels=c("PZAsus-Af2",
                                               "Pza",
                                               "Af2-Unk2",
                                               "Af2",
                                               "Unk2",
                                               "Unk3-Af1",
                                               "Unk3-Unk9",
                                               "Af1-Unk9",
                                               "Unk3",
                                               "Af1",
                                               "Unk9",
                                               "Eu2-Unk4",
                                               "Unk5-Eu2",
                                               "Unk5-Unk4",
                                               "Eu2",
                                               "Unk4",
                                               "Unk5",
                                               "Eu1-Unk7",
                                               "Eu1-Unk6",
                                               "Unk7-Unk6",
                                               "Eu1",
                                               "Unk6",
                                               "Unk7" )) ggplot(table, aes(x = distances, y = Species2)) +
  geom_density_ridges(aes(fill = Species2))
#put un5 and un5 and Eu2 together in one group put un6 and un7 and Eu1 together in another group
un5_4_eu2Gnumbers <-c(un4Gnumbers,un4Gnumbers,eu2Gnumbers) un6_7_eu1Gnumbers <-c(un6Gnumbers,un7Gnumbers,eu1Gnumbers) 
distUn4un5Eu2_Un6un7Eu1<-dist737_pairwisedel[rownames(dist737_pairwisedel) %in% un5_4_eu2Gnumbers, rownames(dist737_pairwisedel) %in% un6_7_eu1Gnumbers] 
matrixDistAll_2=c(distPzaAf2,distPza,distAf2un2,distAf2,distUn2,distAf1un9,distUn3un9,distUn3Af1,distUn3,distAf1,distUn9,distUn4un5Eu2_Un6un7Eu1,distEu2,distUn4,distUn5,distEu1,distUn6,distUn7) 
matrixNameAll_2 <-c(rep("PZAsus-Af2",length(distPzaAf2)),
                  rep("Pza",length(distPza)),
                  rep("Af2-Unk2",length(distAf2un2)),
                  rep("Af2",length(distAf2)),
                  rep("Unk2",length(distUn2)),
                  rep("Af1-Unk9",length(distAf1un9)),
                  rep("Unk3-Unk9",length(distUn3un9)),
                  rep("Unk3-Af1",length(distUn3Af1)),
                  rep("Unk3",length(distUn3)),
                  rep("Af1",length(distAf1)),
                  rep("Unk9",length(distUn9)),
                  rep("Unk4unk5Eu2-Unk6unk7Eu1",length(distUn4un5Eu2_Un6un7Eu1)),
                  rep("Eu2",length(distEu2)),
                  rep("Unk4",length(distUn4)),
                  rep("Unk5",length(distUn5)),
                  rep("Eu1",length(distEu1)),
                  rep("Unk6",length(distUn6)),
                  rep("Unk7",length(distUn7))) table_2 <-data.frame(distances=length(matrixDistAll_2),Species=1,stringsAsFactors = FALSE) for (dist in 1:length(matrixDistAll_2)){
  table_2[dist,1] =matrixDistAll_2[dist]
  table_2[dist,2] =matrixNameAll_2[dist]
}
table_2$Species3 <-factor(table_2$Species,levels=c("PZAsus-Af2",
                                               "Pza",
                                               "Af2-Unk2",
                                               "Af2",
                                               "Unk2",
                                               "Unk3-Af1",
                                               "Unk3-Unk9",
                                               "Af1-Unk9",
                                               "Unk3",
                                               "Af1",
                                               "Unk9",
                                               "Unk4unk5Eu2-Unk6unk7Eu1",
                                               "Eu2",
                                               "Unk4",
                                               "Unk5",
                                               "Eu1",
                                               "Unk6",
                                               "Unk7" )) plot <-ggplot(table_2, aes(x = distances, y = Species3)) +
  geom_density_ridges(aes(fill = Species3))

#Axis of the figures also expressed per kb by diving the number os differences by the size of the alignment used to obtained the genetic distances 
#express the distances in SNP/Kb. The size of the alignment from which snp distances were calculated (alignment_737) was 34308 
250/34308
*1000
[1] 7.28693

500/34308*1000
[1] 14.57386

750/34308*1000
[1] 21.86079
