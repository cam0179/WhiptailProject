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
summary(res)
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

