####  DESeq2 Script for Differential Gene Expression Analysis in 
# whiptail lizards from three different populations
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html



#### Install the DESeq2 package if you do not already have it
## try http:// if https:// URLs are not supported
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

## Load the DESeq2 library 
library(DESeq2)

## Use the Session menu to set working directory to the directory that contains your source files and R scripts
#save your working directory here
setwd("C:/Users/allys/Box/Auburn/FunctionalGenomics/RNAseqPractice/WhiptailProject")

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
### In this pipeline, we are using a gene count matrix that was created by the StringTie mapper.

countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)


### Input the meta data or phenotype data that provides information on where samples came from,
#   what treatment groups they were exposed to, and any other variables you wish to use in analysis. 
coldata <-(read.table("PHENO_DATA.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
#if lines 43 or 45 output "FALSE", your PHENO data file and gene count files rownames are not the same
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
#   Here we want to evaluate differential expression by county samples were collected in
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~county)
#look at it
dds



#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 10 reads mapped.
dds <- dds[ rowSums(counts(dds)) > 10, ]
#check to see how many genes were filtered out
dds

## set factors for statistical analyses
###### Note you need to change condition to whatever treatment you built your model with above


#Lets start with Culberson vs. Brewster
dds$condition <- factor(dds$county, levels=c("Culberson","SantaFe","Brewster"))

################     1.4 Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

# We can order our results table by the smallest adjusted p value:
resOrdered <- res[order(res$padj),]
resOrdered
# We can summarize some basic tallies using the summary function the default is p<0.1.
 
#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


###    1.5.1 MA-plot
## plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
## Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down
plotMA(res, main="DESeq2", ylim=c(-8,8))

#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
# One can then recover the gene identiers by saving the resulting indices:
idx <- identify(res$baseMean, res$log2FoldChange)
# after selecting a gene. You need to press escape to move on
rownames(res)[idx]


##  1.5.2 Plot counts - sanity check!
# You can select the gene to plot by rowname or by numeric index.
plotCounts(dds, gene="gene-IKZF5|IKZF5", intgroup="county")
# You can plot the gene with the lowest adjuated P-value
plotCounts(dds, gene=which.min(res$padj), intgroup="county")
dds

##  Write your results to a file 
write.csv(as.data.frame(resOrdered), file="DGESeq_results_county.csv")  

## 2.1.2 Extracting transformed values
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)

### Heatmap of the count matrix
#library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

library("pheatmap")
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("county", "capture")])
df <- as.data.frame(colData(dds)[,c("county","capture")])
pheatmap(mat, annotation_col = anno)


## 2.2.1 Heatmap of the count matrix
#  library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("county","capture")])
#pheatmap(assay(vsd)[mat,], cluster_rows=TRUE, show_rownames=TRUE,
#        cluster_cols=TRUE, annotation_col=df)
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,  cluster_cols=FALSE, annotation_col=df) 


#2.2.2 Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$county)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# 2.2.3 Principal component plot of the samples
plotPCA(rld, intgroup=c("county"))

plotPCA(rld, intgroup=c("capture"))



##################Using this Analysis results we made these graphs to visualize our resulsts
################## Here we made PCA using ggplot, Volcano graph for our analysis, Venn Diagrams for DGE for each combination and Bargraph to shows up-regulated and down regulated DGEs in each combination

#########################################################################################################################################################################################

#################################### Different types of data visualization methods################################


####get PCA data without making a plot using ggplot

pca<-plotPCA(rld,intgroup="county",ntop=nrow(rld),returnData = TRUE)
ggplot(data=pca,aes_string(x="PC1",y="PC2",color="county"))+geom_point(size=4)+xlab(paste0("PC1:21% variance"))+ylab(paste0("PC2:13%"))+coord_fixed()+scale_color_manual(values=c("#009E73", "#56B4E9", "#E69F00", "#D55E00")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 15, colour = "black"))

####add labels to sample (geom_repel) function
ggplot(data=pca,aes_string(x="PC1",y="PC2",color="county"))+geom_point(size=4)+xlab(paste0("PC1:21% variance"))+ylab(paste0("PC2:13%"))+geom_text_repel(aes(label=name))+coord_fixed()+scale_color_manual(values=c("#009E73", "#56B4E9", "#E69F00", "#D55E00")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 15, colour = "black"))
##################################################################################################################################################################################

#######Volcano plot for differential expression gene data

#install enhanced volcano package and textshaping package
BiocManager::install("EnhancedVolcano")
install.packages("textshaping")

#Open the Enhancedvolcano library
library(EnhancedVolcano)
library(textshaping)
write.csv(as.data.frame(res),file = "volcano.csv")

#Manually put a name to first column in volcano.csv

#open the volcano csv file as a dataframe
res_volcano<-as.data.frame(read.csv("volcano.csv"))
res_volcano

#use EnhancedVolcano to draw the volcano plot 
EnhancedVolcano(res_volcano,x = "log2FoldChange", y = "padj",lab = res_volcano$Gene_Name ,labSize = 2)

################################################################################################################################################

########Venn diagram for differently Differently expressed genes

########Lets start with SantaFe vs. Brewster
dds$condition <- factor(dds$county, levels=c("Culberson","SantaFe","Brewster"))
dds <- DESeq(dds)
res_SF_bre <- results(dds)
res_SF_bre

#oredered the genes according to Padj values
res_SF_bre_Ordered <- res_SF_bre[order(res_SF_bre$padj),]
res_SF_bre_Ordered

#only extract the genes with lower 0.05 padj in SantaFe vs. Brewster
res_SF_bre_Ordered_0.05 <- subset(res_SF_bre_Ordered, padj < 0.05)

#make a excel file for this data
write.csv(as.data.frame(res_SF_bre_Ordered_0.05), file="DGE_SantaFE_vs_Brewster.csv")

###########Next Culberson vs Brewster

dds$condition <- factor(dds$county, levels=c("SantaFe","Culberson","Brewster"))

#here we relevel this to get culberson vs brewster DEGs
dds$county<-relevel(dds$county,ref = "Brewster")
dds <- DESeq(dds)
res_cul_bre <- results(dds)
res_cul_bre

#oredered the genes according to Padj values
res_cul_bre_Ordered <- res_cul_bre[order(res_cul_bre$padj),]
res_cul_bre_Ordered

#only extract the genes with lower 0.05 padj in SantaFe vs. Brewster
res_cul_bre_Ordered_0.05 <- subset(res_cul_bre_Ordered, padj < 0.05)

#make a excel file for this data
write.csv(as.data.frame(res_cul_bre_Ordered_0.05), file="DGE_Culberson_vs_Brewster.csv")


###########Next Culberson vs Santafe

dds$condition <- factor(dds$county, levels=c("SantaFe","Culberson","Brewster"))

#here we relevel this to get culberson vs brewster DEGs
dds$county<-relevel(dds$county,ref = "SantaFe")
dds <- DESeq(dds)
res_cul_San <- results(dds)
res_cul_San

#oredered the genes according to Padj values
res_cul_San_Ordered <- res_cul_San[order(res_cul_San$padj),]
res_cul_San_Ordered

#only extract the genes with lower 0.05 padj in SantaFe vs. Brewster
res_cul_San_Ordered_0.05 <- subset(res_cul_San_Ordered, padj < 0.05)

#make a excel file for this data
write.csv(as.data.frame(res_cul_San_Ordered_0.05), file="DGE_Culberson_vs_SantaFe.csv")


####now we have 3 excel files that include DEG for "culberson vs SantaFe",Culberson vs Brester", and "SantaFe vs Brewster".
####Using DGE names in these file we can make a venn diagram that shows differences and similarities in these situations.
####To make this ven diagram we used Venny 2.1.0 online application. https://bioinfogp.cnb.csic.es/tools/venny/

#########################################################################################################################################

#######################Making a bar graph that showing numbers of up regulated and down regulated genes

#First check Santa Fe vs Brewster
SanVSBre<-read.csv("DGE_SantaFE_vs_Brewster.csv")

#get the number of upregulated genes in SantaFe vs Brewster
SanVSBreUP<-subset(SanVSBre, log2FoldChange >0)
SanVSBreUP
nrow(SanVSBreUP)

#get the number of downegulated genes in SantaFe vs Brewster
SanVSBreDown<-subset(SanVSBre, log2FoldChange <0)
SanVSBreDown
nrow(SanVSBreDown)

##########################

####Next check Culberson vs Brewster
CulVSBre<-read.csv("DGE_Culberson_vs_Brewster.csv")

#get the number of upregulated genes in SantaFe vs Brewster
CulVSBreUP<-subset(CulVSBre, log2FoldChange >0)
CulVSBreUP
nrow(CulVSBreUP)

#get the number of downegulated genes in SantaFe vs Brewster
CulVSBreDown<-subset(CulVSBre, log2FoldChange <0)
CulVSBreDown
nrow(CulVSBreDown)

#############Next Culberson vs Santa Fe

####Next check Culberson vs Brewster
CulVSSF<-read.csv("DGE_Culberson_vs_SantaFe.csv")

#get the number of upregulated genes in SantaFe vs Brewster
CulVSSFUP<-subset(CulVSSF, log2FoldChange >0)
CulVSSFUP
nrow(CulVSSFUP)

#get the number of downegulated genes in SantaFe vs Brewster
CulVSSFDown<-subset(CulVSSF, log2FoldChange <0)
CulVSSFDown
nrow(CulVSSFDown)

########Here all the nrow() values recorded in a table in Excel. Finally using Excel we made the bargraph

