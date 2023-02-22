## ----Install edgeR, message=FALSE------------------------------------------------------------------------------------------------------------------------------------
#BiocManager::install("edgeR")


## ----Select edgeR----------------------------------------------------------------------------------------------------------------------------------------------------
library(edgeR)


## ----Data Load-------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("~/Documents/BIOTECHIE/Unimi/Methods in Bioinformatics/info18/files_method_bioinformatics_exam")
raw_counts <- read.csv("data.csv", header = T,  sep = "\t", row.names = 1) 
group <- read.csv("conditions.csv", header = T, sep = "\t", row.names = 1)
#row.names=1 because the objects of the first column give the names to the different rows of the data frame
#header= TRUE because the objects of the first row give the names to the different columns of the data frame


## ----Order-----------------------------------------------------------------------------------------------------------------------------------------------------------
#with this line of code we create a dataframe which will have as columns names the names of the rows in the data.frame 'group' in the SAME order
counts <- raw_counts[,rownames(group)] 
#this instead is a logical vector that demonstrate that what we have done is correct
colnames(counts) == rownames(group)


## ----DGEList---------------------------------------------------------------------------------------------------------------------------------------------------------
#we set that the library size of each replicate (counts columns) is equal to the sum of each value whithin the column i.e. the reads of each gene for each condition
countsDGE <- DGEList(counts= counts, group = group$tissue, remove.zeros = T, lib.size =colSums(counts))


## ----Filter----------------------------------------------------------------------------------------------------------------------------------------------------------
#we have summarized this operation in one string of code
countsDGE$counts <- countsDGE$counts[apply(countsDGE$counts>10,1,sum)>0,]


## ----TMM-------------------------------------------------------------------------------------------------------------------------------------------------------------
#this will add a column named "norm.factors" ro the samples data.frame of the DEGList
N.countsDGE <- calcNormFactors(countsDGE, method = "TMM")


## ----MDS-------------------------------------------------------------------------------------------------------------------------------------------------------------

DiffCol <- c("red","orange","blue","green","mediumorchid1","black")[N.countsDGE$samples$group]
plotMDS(N.countsDGE, col=DiffCol, pch = 19)
legend("topright",fill=c("red","orange","blue","green","mediumorchid1","black"),legend=levels(N.countsDGE$samples$group))



## ----NB Fitting------------------------------------------------------------------------------------------------------------------------------------------------------
#if not specified the type of dispersion used will be tagwise
D.countsDGE <- estimateDisp(N.countsDGE)
#by using exactTest an object DGEExact is returned
eT <- exactTest(D.countsDGE, pair = c("Frontal Cortex (BA9)","Substantia nigra"))



## ----not Filtered DEGs-----------------------------------------------------------------------------------------------------------------------------------------------
Ctt <- topTags(eT, n=nrow(eT)) 
#these are the genes that are presumably DE but that are not yet filtered for an FDR


## ----Filtered Degs---------------------------------------------------------------------------------------------------------------------------------------------------
#real DEGs have an FDR lower than 0.01
DE_UP <- (Ctt[Ctt$table$FDR <= 0.01 & Ctt$table$logFC>0,])
DE_DOWN <- (Ctt[Ctt$table$FDR <= 0.01 & Ctt$table$logFC<0,])
#genes that have an FDR higher than 0.01 are not considered DEGs
notDE_UP <- (Ctt[Ctt$table$FDR > 0.01 & Ctt$table$logFC>0,])
notDE_DOWN <- (Ctt[Ctt$table$FDR > 0.01 & Ctt$table$logFC<0,])


## ----DEGs Boxplots---------------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))#to visualize the boxplots in a 2x2 grid
boxplot(DE_UP$table$logFC, main = "logFC DE_UP", outline = F, col = 'mediumorchid')
boxplot(DE_DOWN$table$logFC, main = "logFC DE_DOWN", outline = F, col = 'mediumorchid2')
boxplot(notDE_UP$table$logFC, main = "logFC notDE_UP", outline = F, col = 'mediumorchid4')
boxplot(notDE_DOWN$table$logFC, main = "logFC notDE_DOWN", outline = F, col = 'violet')


## ----S.Nigra Fr. Cortex filter---------------------------------------------------------------------------------------------------------------------------------------
#first we filter the conditions of the replicates for which neither Frontal Cortex nor Substantia Nigra is sampled
group2 <- group[-c(1:12,19:30),] 
#then we filter out gene counts that are referring only to the replicates for which our tissues are sampled by indexing the columns of the columns of the counts data.frame
counts2 <- counts[,rownames(group2)]


## ----Subsetting Fr. Cortex-------------------------------------------------------------------------------------------------------------------------------------------
#Subsetting Frontal Cortex
groupFrC <- group2[-c(7:12),]
countsFrC <- counts2[,rownames(groupFrC)]



## ----Filtering-------------------------------------------------------------------------------------------------------------------------------------------------------
#we need to set group = groupFrC$sex to add to the DEGList object a column in the samples data.frame that assigns the sex membership of each sample
DGEFrC <- DGEList(counts = countsFrC, group = groupFrC$sex)
DGEFrC$counts <- DGEFrC$counts[apply(DGEFrC$counts>10,1,sum)>1,]


## ----Normalization---------------------------------------------------------------------------------------------------------------------------------------------------
#normalization
DGEFrC <- calcNormFactors(DGEFrC)
#pseudocounts are essential to deal with 0s
logDGEFrC <- log(DGEFrC$counts + 0.01)


## ----PCA-------------------------------------------------------------------------------------------------------------------------------------------------------------
#PCA 
col.sex <- c("pink","lightblue")[DGEFrC$samples$group]
#it is essential to apply the prcomp() function on the columns of the data.frame containing the counts because the columns contain the id of each replicate; i need to compare the replicates, not the genes
PCAFrC <- prcomp(t(logDGEFrC))
plot(PCAFrC$x, col = col.sex, pch = 19, main = "Frontal Cortex PCA")
legend("bottomright",fill = c("pink","lightblue"), legend = levels(DGEFrC$samples$group))


## ----NB Fitting Fr. Cortex-------------------------------------------------------------------------------------------------------------------------------------------
D.DGEFrC <- estimateDisp(DGEFrC)


## ----DEGs------------------------------------------------------------------------------------------------------------------------------------------------------------
eTFrC <- exactTest(D.DGEFrC, pair = c(1,2))


## ----filtered DEGs---------------------------------------------------------------------------------------------------------------------------------------------------
FrCTT <- topTags(eTFrC, n=nrow(eTFrC))
#in task 7 we filtered data both for logFC and for FDR, here we filter only for an FRD
genes_sex_spec <- FrCTT[FrCTT$table$FDR <= 0.05,]
#header genes
head(genes_sex_spec)


## ----used genes legend-----------------------------------------------------------------------------------------------------------------------------------------------
library(biomaRt)
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes_used <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id",values = rownames(genes_sex_spec), mart = ensembl)
genes_used


## ----all about Substantia Nigra--------------------------------------------------------------------------------------------------------------------------------------

#Subsetting Substantia Nigra

groupSubNig <- group2[c(7:12),]
countsSubNig <- counts2[,rownames(groupSubNig)]

#DGEList + filtering
DGESubNig <- DGEList(counts = countsSubNig, group = groupSubNig$sex)
DGESubNig$counts <- DGESubNig$counts[apply(DGESubNig$counts>10,1,sum)>1,]

#Normalization TMM + log scaling
DGESubNig <- calcNormFactors(DGESubNig)
logDGESubNig <- log(DGESubNig$counts + 0.01)

#PCA

col.sex <- c("pink","lightblue")[DGESubNig$samples$group]
PCASubNig <- prcomp(t(logDGESubNig))

plot(PCASubNig$x, col = col.sex, pch = 19, main = "Substantia Nigra's PCA")
legend("bottomright", fill = c("pink","lightblue"), legend = levels(DGESubNig$samples$group))

#Diff expression analysis to extract sex specific genes in Substantia Nigra
#as in part 1 we need at first to estimate dispersion. We already computed normalization and applied filtering as requested.

D.DGESubNig <- estimateDisp(DGESubNig)

#Fisher exact test

eTSubNig <- exactTest(D.DGESubNig, pair = c(1,2))

#After Fisher exact test, we can proceed with considering all genes that show FDR <= 0.05
#as sex specific DEGs

#Let's compute TopTags

SubNigTT <- topTags(eTSubNig, n=nrow(eTSubNig))

#... and then filter out our genes

genes_sex_spec2 <- SubNigTT[SubNigTT$table$FDR <= 0.05,]

#header genes
head(genes_sex_spec2)
genes_used2 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id",values = rownames(genes_sex_spec2), mart = ensembl)
genes_used2


## ----Venn for common sex---------------------------------------------------------------------------------------------------------------------------------------------
library('VennDiagram')
Frc_string <- rownames(genes_sex_spec)
SubNig_string <- rownames(genes_sex_spec2)
listalo <- list(Frc_string, SubNig_string)
venn.diagram(listalo, category.names = c("F.Cortex", "S.Nigra"),filename = "DEGsex.png", imagetype = "png",main = "Sex specific genes in Frontal Cortex and Substantia Nigra", col = c("mediumorchid", "mediumorchid4"))


## ----Venn DEGs and sex DEGs------------------------------------------------------------------------------------------------------------------------------------------
#VENN DIAGRAM between DEGs in common part (common part), S.Nigra and Fr. Cortex
Ctt_string <- rownames(Ctt[Ctt$table$FDR <= 0.01,])
Camillibus <- list(Frc_string, SubNig_string, Ctt_string)

venn.diagram(Camillibus, category.names = c("F.Cortex", "S.Nigra", "Common part"),filename = "j@mpaolo.png", imagetype = "png", main = "All DEGs vs Sex Specific DEGs of S.Nigra and F.Cortex",col = c("mediumorchid", "mediumorchid4", "violet"))


