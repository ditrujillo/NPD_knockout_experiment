##### Clean the working space
rm (list = ls ())

##### (Install packages) --------------------------------------------------
source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq")
a
biocLite("GOstats")
biocLite("edgeR")

##### Load packages -------------------------------------------------------
library(DESeq)
library(GOstats)
library(GO.db)
if (!require("gplots")) { install.packages("gplots", dependencies = TRUE) }
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) }
if (!require("reshape2")) { install.packages("reshape2", dependencies = TRUE) }
if (!require("edgeR")) { install.packages("edgeR", dependencies = TRUE) }
if (!require("limma")) { install.packages("limma", dependencies = TRUE) }

##### Import Data ---------------------------------------------------------
#Groups: HM340=WT, Gm12-7=NPD24, Gm8-5=NPD245, Gm26-1=NPD124
sampleFiles <- c ("HM340_N1_S40_counts.txt", "HM340_N3_S44_counts.txt", "HM340_N4_S48_counts.txt","Gm12-7_N1_S41_counts.txt", "Gm12-7_N3_S45_counts.txt", "Gm12-7_N4_S49_counts.txt", "Gm8-5_N1_S43_counts.txt", "Gm8-5_N3_S47_counts.txt", "Gm8-5_N4_S51_counts.txt", "Gm26-1_N1_S42_counts.txt", "Gm26-1_N3_S46_counts.txt", "Gm26-1_N4_S50_counts.txt")
sampleName <- c ("WT_1", "WT_2", "WT_3", "NPD24_1", "NPD24_2", "NPD24_3", "NPD245_1", "NPD245_2", "NPD245_3", "NPD124_1", "NPD124_2", "NPD124_3")
group <- factor(c(rep("WT",3), rep("NPD24",3), rep("NPD245",3), rep("NPD124",3)))
sampleTable <- data.frame(sampleName = sampleName, fileName = sampleFiles, group=group)
directory <- c("~/Data/PDPs/HTseq/")
cds <- newCountDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory)

##### Normalize the data --------------------------------------------------
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
norm<-counts(cds, normalized=TRUE)

##### Test for DE genes ---------------------------------------------------
testNPD24 = nbinomTest(cds, "WT", "NPD24")
testNPD245 = nbinomTest(cds, "WT", "NPD245")
testNPD124 = nbinomTest(cds, "WT", "NPD124")
testNPD24.NPD245 = nbinomTest(cds, "NPD24", "NPD245")
testNPD24.NPD124 = nbinomTest(cds, "NPD24", "NPD124")
testNPD124.NPD245 = nbinomTest(cds, "NPD124", "NPD245")

##### (Look at specific genes) --------------------------------------------
subset(testNPD24, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))
#MtCAM1=Medtr3g055570, MtCAM2=Medtr3g055585, MtCAM3=Mtr.37968.1.S1_at, MtCAM5=Medtr3g055510
subset(testNPD24, id=="Medtr3g067615")
subset(testNPD245, id=="Medtr3g067615")
subset(testNPD124, id=="Medtr3g067615")
subset(testNPD24.NPD245, id=="Medtr3g106480")
subset(testNPD24.NPD124, id=="Medtr3g106480")
subset(testNPD124.NPD245, id=="Medtr7g109920")

##### (Plot correlations of normalized expression levels) -----------------
cor<-cor(norm, use="pairwise.complete.obs", method="pearson")
#pdf("heatmap.pdf", pointsize=10, font="Helvetica")
heatmap.2(cor, trace="none", col="redgreen", density.info="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", margin=c(6,6))
#dev.off()

hist(cor, breaks=50)
dist<-as.dist((1-cor)/2)
#pdf(file="tree.pdf")
plot(hclust(dist))
#dev.off()

##### (Look at the number and plot DE genes) ------------------------------
nrow(subset(test, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2)))
test=testNPD24
#Plot genes with p-val<0.05 and p-adj<0.05
plot(test$log2FoldChange~log10(test$baseMean), xlab="log10 Mean Gene Expression", ylab = "log2 Fold Change", main="HM340 vs. NPD1/2/4")
points(test$log2FoldChange[test$pval<0.05]~log10(test$baseMean[test$pval<0.05]), col="red", pch=16)
points(test$log2FoldChange[test$padj<0.05]~log10(test$baseMean[test$padj<0.05]), col="blue", pch=16)     

##### Choose sets of DE genes with p-adjusted<0.1 ----------------------------------------
#Padj < 0.1, greater than 3-fold change
Sub.Ad24 = subset(testNPD24, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Ad245 = subset(testNPD245, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Ad124 = subset(testNPD124, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Ad24.245 = subset(testNPD24.NPD245, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Ad24.124 = subset(testNPD24.NPD124, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Ad124.245 = subset(testNPD124.NPD245, padj<0.1 & !(foldChange >= 1/3 & foldChange <= 3))$id
#Create list of all DE genes. Remove targeted NPD genes, which show up as DE
AllDE = sort(Reduce(union, list(Sub.Ad24,Sub.Ad245,Sub.Ad124,Sub.Ad24.245,Sub.Ad24.124,Sub.Ad124.245)))
AllDE = AllDE[AllDE!="Medtr2g103303" & AllDE!="Medtr2g103307" & AllDE!="Medtr2g103330"]
length(AllDE)
#Write files  with All DE genes, and specific tests
setwd("~/Data/PDPs/DESeq")
write(AllDE, file = "AllDE.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Ad24, file = "Sub.Ad24.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Ad245, file = "Sub.Ad245.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Ad124, file = "Sub.Ad124.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Ad24.245, file = "Sub.Ad24.245.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Ad24.124, file = "Sub.Ad24.124.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Ad124.245, file = "Sub.Ad124.245.list", ncolumns = 1, append = FALSE, sep = "\n")

##### Choose sets of DE genes with p-val<0.01 --------------------------------
Sub.Pv24 = subset(testNPD24, pval<0.01 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Pv245 = subset(testNPD245, pval<0.01 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Pv124 = subset(testNPD124, pval<0.01 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Pv24.245 = subset(testNPD24.NPD245, pval<0.01 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Pv24.124 = subset(testNPD24.NPD124, pval<0.01 & !(foldChange >= 1/3 & foldChange <= 3))$id
Sub.Pv124.245 = subset(testNPD24.NPD124, pval<0.01 & !(foldChange >= 1/3 & foldChange <= 3))$id
#Create list of all DE genes. Remove targeted NPD genes, which show up as DE
AllDEPv = sort(Reduce(union, list(Sub.Pv24,Sub.Pv245,Sub.Pv124,Sub.Pv24.245,Sub.Pv24.124,Sub.Pv124.245)))
AllDEPv = AllDEPv[AllDEPv!="Medtr2g103303" & AllDEPv!="Medtr2g103307" & AllDEPv!="Medtr2g103330"]
length(AllDEPv)
#Write files  with All DE genes, and specific tests
setwd("~/Data/PDPs/DESeq")
write(AllDEPv, file = "AllDEPv.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv24, file = "Sub.Pv24.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv245, file = "Sub.Pv245.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv124, file = "Sub.Pv124.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv24.245, file = "Sub.Pv24.245.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv24.124, file = "Sub.Pv24.124.list", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv124.245, file = "Sub.Pv124.245.list", ncolumns = 1, append = FALSE, sep = "\n")

##### Load sets of DE genes ----------------------------------------------
setwd("~/Data/PDPs/DESeq")
AllDE = scan("AllDE.list", what="", sep="\n")
Sub.Ad24 = scan("Sub.Ad24.list", what="", sep="\n")
Sub.Ad245 = scan("Sub.Ad245.list", what="", sep="\n")
Sub.Ad124 = scan("Sub.Ad124.list", what="", sep="\n")
Sub.Ad24.245 = scan("Sub.Ad24.245.list", what="", sep="\n")
Sub.Ad24.124 = scan("Sub.Ad24.124.list", what="", sep="\n")
Sub.Ad124.245 = scan("Sub.Ad124.245.list", what="", sep="\n")
AllDEPv = scan("~/Data/PDPs/DESeq/AllDEPv.list", what="", sep="\n")
Sub.Pv24 = scan("Sub.Pv24.list", what="", sep="\n")
Sub.Pv124 = scan("Sub.Pv124.list", what="", sep="\n")
Sub.Pv245 = scan("Sub.Pv245.list", what="", sep="\n")
Sub.Pv24.124 = scan("Sub.Pv24.124.list", what="", sep="\n")
Sub.Pv24.245 = scan("Sub.Pv24.245.list", what="", sep="\n")
Sub.Pv124.245 = scan("Sub.Pv24.245.list", what="", sep="\n")
data.frame(norm[row.names(norm) %in% Sub.Ad24, ])

##### Heatmap of ALL DE genes with padj<0.1 -------------------------------
DEtable = data.frame(norm[row.names(norm) %in% AllDE, ])
ScaleDE = t(scale(t(DEtable)))
#Yellow indicates high  expression while red indicates low expression
heatmap.2(ScaleDE, dendrogram='row', Rowv=TRUE, Colv=FALSE, margin=c(6,6), 
          key=FALSE, trace='none',density.info="none",lhei=c(0.1,8), lwid=c(1,1))
venn(list(Sub.Ad24,Sub.Ad245,Sub.Ad124,Sub.Ad24.245,Sub.Ad24.124))

########## Figure 5 ########## 
##### Heatmap of DE genes (padj<0.1) compared to WT or NPD24 --------------
### Choose subsets of DE genes and assign colors --------------------------
## A gene appears as DE when compared against at least HM340(WT)
## or NPD2/4, which had a WT nodule phenotype
# Implication: Even though the adjusted p-value cutoff is 0.1, 
# if a gene appears as DE in multiple comparisons, this gives 
# us more confidence that it truly is DE

## Genes that are DE with padj<0.1 compared to WT or NPD24 (117 genes)
DEvsWTorNPD24=sort(Reduce(union, list(Sub.Ad24,Sub.Ad245,Sub.Ad124,Sub.Ad24.245,Sub.Ad24.124)))
DEvsWTorNPD24 = DEvsWTorNPD24[DEvsWTorNPD24!="Medtr2g103303" & DEvsWTorNPD24!="Medtr2g103307" & DEvsWTorNPD24!="Medtr2g103330"]
length(DEvsWTorNPD24)
write(DEvsWTorNPD24, file = "DEvsWTorNPD24.list", ncolumns = 1, append = FALSE, sep = "\n")
## Genes that are DE in NPD2/4, in comparison to WT (p-adj<0.1) 
CommonDE24 = Sub.Ad24
## For NPD2/4/5 and NPD1/2/4, keep only the ones that appear in multiple 
## comparisons (p-val<0.01)
## Genes that are DE in NPD2/4/5, in comparison to both WT or NPD2/4 (p-val<0.01) 
CommonDE245 = intersect(Sub.Pv245, Sub.Pv24.245)
## Genes that are DE in NPD1/2/4, in comparison to both WT or NPD2/4 (p-val<0.01) 
CommonDE124 = intersect(Sub.Pv124, Sub.Pv24.124)
## A gene is DE in NPD1/2/4 and NPD2/4/5 relative to WT and NPD2/4:
CommonDE124.245 = intersect(CommonDE245,CommonDE124)
## Final list of genes with p-adj<0.1
CommonDEs = intersect(DEvsWTorNPD24,sort(Reduce(union, list(CommonDE24,CommonDE245,CommonDE124,CommonDE124.245))))
write(CommonDEs, file = "DEmin2vsWTorNPD24.list", ncolumns = 1, append = FALSE, sep = "\n")
length(CommonDEs)
DEtable = data.frame(norm[row.names(norm) %in% CommonDEs, ])
ScaleDE = t(scale(t(DEtable)))
## Add gene categories to table
DEtable$Categories <- "Cat0"
DEtable$Categories[row.names(DEtable) %in% CommonDE245] <- "Cat1"
DEtable$Categories[row.names(DEtable) %in% CommonDE124] <- "Cat2"
DEtable$Categories[row.names(DEtable) %in% CommonDE124.245] <- "Cat4"
DEtable$Categories[row.names(DEtable) %in% CommonDE24] <- "Cat3"
length(DEtable$Categories=="Cat2")
## Color DE genes by categories
# Assign colors: Red = CommonDE245; 
#                Green = CommonDE124
#                Blue = CommonDE24
#                Purple = CommonDE124.245
colours <- ifelse(DEtable$Categories == "Cat1", "red", 
                  ifelse(DEtable$Categories == "Cat2", "darkgreen", 
                         ifelse(DEtable$Categories == "Cat3", "blue", 
                                ifelse(DEtable$Categories == "Cat4", "purple",
                                       "white"))))
### Create heatmap --------------------------------------------------------
png(filename='~/Data/PDPs/DESeq/DEvsWTorNPD24.png', width=1200, height=1000)
heatmap.2(ScaleDE, dendrogram='row', Rowv=TRUE, Colv=FALSE, margin=c(11,83), 
          key=FALSE, trace='none',density.info="none",lhei=c(0.1,10), lwid=c(1,6),
          cexRow = 1.8, cexCol = 2, RowSideColors=colours)
graphics.off()
### Create heatmap with annotations ---------------------------------------
# Get heatmap order
hm = heatmap.2(ScaleDE,dendrogram='row', Rowv=TRUE, Colv=FALSE)
hclust = as.hclust(hm$rowDendrogram)
geneorder = cutree(hclust, h=10)[hclust$order]
# Get genes from annotation file
Annotations = read.delim("~/Data/Mt4.0/Data/genome/MtDescriptionList.txt", header=FALSE)
index = match(labels(geneorder), Annotations$V1)
annotationorder = Annotations$V2[rev(index)]
GeneTable = subset(Annotations, V1 %in% labels(geneorder))
GeneTable = GeneTable[rev(match(labels(geneorder), GeneTable$V1)),]
# Add labels to heatmap
y <- length(labels(geneorder)):1
#x <- 1:length(labels(geneorder))
my_plot1 <- function(data){
  heatmap.2(data, dendrogram='row', Rowv=TRUE, Colv=FALSE, margin=c(11,105),
            key=FALSE, trace='none',density.info="none",lhei=c(0.8,10), lwid=c(0.5,10),
            cexRow = 1.5, cexCol = 1.5, RowSideColors=colours)
  text(x=0.22, y=0.028+0.948*y/length(y), labels=annotationorder, cex=1, pos=4, col="darkblue")
}
png(filename='~/Data/PDPs/DESeq/DEvsWTorNPD24Annotated.png', width=1200, height=1000)
my_plot1(ScaleDE)
graphics.off()

##### Heatmap of DE genes (padj<0.1) in at least two samples --------------
### Choose subsets of DE genes and assign colors --------------------------
## A gene appears as DE when compared against at least 2 lines:
# Implication: Even though the adjusted p-value cutoff is 0.1, 
# if a gene appears as DE in multiple comparisons, this gives 
# us more confidence that it truly is DE
DEmin2=sort(Reduce(union, list(intersect(Sub.Ad24, Sub.Ad245),
                               intersect(Sub.Ad24, Sub.Ad124),
                               intersect(Sub.Ad24, Sub.Ad24.245),
                               intersect(Sub.Ad24, Sub.Ad24.124),
                               intersect(Sub.Ad24, Sub.Ad124.245),
                               intersect(Sub.Ad245, Sub.Ad124),
                               intersect(Sub.Ad245, Sub.Ad24.245),
                               intersect(Sub.Ad245, Sub.Ad24.124),
                               intersect(Sub.Ad245, Sub.Ad124.245),
                               intersect(Sub.Ad124, Sub.Ad24.245),
                               intersect(Sub.Ad124, Sub.Ad24.124),
                               intersect(Sub.Ad124, Sub.Ad124.245),
                               intersect(Sub.Ad24.245, Sub.Ad24.124),
                               intersect(Sub.Ad24.245, Sub.Ad124.245),
                               intersect(Sub.Ad24.124, Sub.Ad124.245))))
DEmin2 = DEmin2[DEmin2!="Medtr2g103303" & DEmin2!="Medtr2g103307" & DEmin2!="Medtr2g103330"]
length(DEmin2)
write(DEmin2, file = "DEmin2.list", ncolumns = 1, append = FALSE, sep = "\n")
DEtable = data.frame(norm[row.names(norm) %in% DEmin2, ])
ScaleDE = t(scale(t(DEtable)))
## Genes that are DE in NPD2/4/5, in comparison to at least two other lines: 
# WTvs.NPD2/4/5, NPD2/4vs.NPD2/4/5, NPD1/2/4vs.NPD2/4/5
CommonDE245=sort(Reduce(union, list(intersect(Sub.Ad245, Sub.Ad24.245),
                                    intersect(Sub.Ad245, Sub.Ad124.245),
                                    intersect(Sub.Ad24.245, Sub.Ad124.245))))
DEtable$Categories <- "Cat0"
DEtable$Categories[row.names(DEtable) %in% CommonDE245] <- "Cat1"
## Genes that are DE in NPD1/2/4, in comparison to at least two other lines: 
# WTvs.NPD1/2/4, NPD2/4vs.NPD1/2/4, NPD1/2/4vs.NPD2/4/5
CommonDE124=sort(Reduce(union, list(intersect(Sub.Ad124, Sub.Ad24.124),
                                    intersect(Sub.Ad124, Sub.Ad124.245),
                                    intersect(Sub.Ad24.124, Sub.Ad124.245))))
DEtable$Categories[row.names(DEtable) %in% CommonDE124] <- "Cat2"
## Genes that are DE in NPD2/4, in comparison to at least two other lines: 
# WTvs.NPD2/4, NPD2/4vs.NPD1/2/4, NPD2/4vs.NPD2/4/5, plus Medtr3g106480, which was close to the cutoff
CommonDE24=sort(Reduce(union, list(intersect(Sub.Ad24, Sub.Ad24.124),
                                   intersect(Sub.Ad24, Sub.Ad24.245),
                                   intersect(Sub.Ad24.124, Sub.Ad24.245), "Medtr3g106480")))
DEtable$Categories[row.names(DEtable) %in% CommonDE24] <- "Cat3"
## A gene is DE in NPD1/2/4 and NPD2/4/5 relative to WT and NPD2/4:
CommonDE124.245= intersect(intersect(Sub.Ad124, Sub.Ad24.124),
                           intersect(Sub.Ad245, Sub.Ad24.245))
DEtable$Categories[row.names(DEtable) %in% CommonDE124.245] <- "Cat4"
## Color DE genes by categories
# Assign colors: Red = CommonDE245; 
#                Green = CommonDE124
#                Blue = CommonDE24
#                Purple = CommonDE124.245
colours <- ifelse(DEtable$Categories == "Cat1", "red", 
                  ifelse(DEtable$Categories == "Cat2", "darkgreen", 
                         ifelse(DEtable$Categories == "Cat3", "blue", 
                                ifelse(DEtable$Categories == "Cat4", "purple",
                                       "white"))))
### Create heatmap --------------------------------------------------------
png(filename='~/Data/PDPs/DESeq/DEmin2.png', width=1200, height=1000)
heatmap.2(ScaleDE, dendrogram='row', Rowv=TRUE, Colv=FALSE, margin=c(11,83), 
          key=FALSE, trace='none',density.info="none",lhei=c(0.1,10), lwid=c(1,6),
          cexRow = 1.8, cexCol = 2, RowSideColors=colours)
graphics.off()
### Create heatmap with annotations ---------------------------------------
# Get Heatmap order
hm = heatmap.2(ScaleDE,dendrogram='row', Rowv=TRUE, Colv=FALSE)
hclust = as.hclust(hm$rowDendrogram)
geneorder = cutree(hclust, h=10)[hclust$order]
# Get genes from annotation file
Annotations = read.delim("~/Data/Mt4.0/Data/genome/MtDescriptionList.txt", header=FALSE)
index = match(labels(geneorder), Annotations$V1)
annotationorder = Annotations$V2[rev(index)]
GeneTable = subset(Annotations, V1 %in% labels(geneorder))
# Add labels to heatmap
y <- length(labels(geneorder)):1
#x <- 1:length(labels(geneorder))
my_plot1 <- function(data){
  heatmap.2(data, dendrogram='row', Rowv=TRUE, Colv=FALSE, margin=c(11,105),
            key=FALSE, trace='none',density.info="none",lhei=c(0.8,10), lwid=c(0.5,10),
            cexRow = 1.5, cexCol = 1.5, RowSideColors=colours)
  text(x=0.22, y=0.028+0.948*y/length(y), labels=annotationorder, cex=1, pos=4, col="darkblue")
}

png(filename='~/Data/PDPs/DESeq/DEmin2Annotated.png', width=1200, height=1000)
my_plot1(ScaleDE)
graphics.off()

##### Gene Ontology Enrichment  ---------------------------------------------
### Load gene list to test ------------------------------------------------
#Get genes to test
setwd("~/Data/PDPs/DESeq")
#Output: AllDE.enrichment.txt, NPD24.enrichment.txt, NPD245.enrichment.txt, NPD124.enrichment.txt
#TestGenes = scan("AllDE.list", what="", sep="\n")
#TestGenes = scan("Sub.Ad24.list", what="", sep="\n")
#TestGenes = sort(Reduce(union, list(scan("Sub.Ad245.list", what="", sep="\n"),scan("Sub.Ad24.245.list", what="", sep="\n"))))
#TestGenes = sort(Reduce(union, list(scan("Sub.Ad124.list", what="", sep="\n"),scan("Sub.Ad24.124.list", what="", sep="\n"))))

#Output: AllDEPv.enrichment.txt, NPD24Pv.enrichment.txt, NPD245Pv.enrichment.txt, NPD124Pv.enrichment.txt
#TestGenes = scan("~/Data/PDPs/DESeq/AllDEPv.list", what="", sep="\n")
#TestGenes = scan("Sub.Pv24.list", what="", sep="\n")
#TestGenes = sort(Reduce(union, list(scan("Sub.Pv124.list", what="", sep="\n"),scan("Sub.Pv24.124.list", what="", sep="\n"))))
#TestGenes = sort(Reduce(union, list(scan("Sub.Pv245.list", what="", sep="\n"),scan("Sub.Pv24.245.list", what="", sep="\n"))))

#Output: AllDEPv.enrichment.txt, NPD24Pv.enrichment.txt, NPD245Pv.enrichment.txt, NPD124Pv.enrichment.txt
#TestGenes = intersect(scan("Sub.Pv124.list", what="", sep="\n"),scan("Sub.Pv24.124.list", what="", sep="\n"))
#TestGenes = intersect(scan("Sub.Pv245.list", what="", sep="\n"),scan("Sub.Pv24.245.list", what="", sep="\n"))

#Output: DEvsWTorNPD24.enrichment.txt
TestGenes = scan("~/Data/PDPs/DESeq/DEvsWTorNPD24.list", what="", sep="\n")
### GO enrichment Tests ---------------------------------------------------
setwd("~/Data/PDPs/GO")
source('~/Data/PDPs/GO/perform_enrichment.R')
#Matrix backbone
AllMt = scan("~/Data/Mt4.0/Data/genome/UniqueList.txt", what="", sep="\n")
GeneMAT = matrix(data = 0, nrow = length(AllMt), ncol = 1)
row.names(GeneMAT) = c(AllMt)
colnames(GeneMAT) = c("DE")
#Put DE gene list into a matrix format. 1 is differentially expressed
for (i in 1:length(TestGenes)){
  temp <- which(TestGenes[i] == row.names(GeneMAT))
  GeneMAT[temp,1] <- 1}
##Medicago GO term propagation
#GOF is Medicago go term propagation file
GOF <- read.delim("~/Data/PDPs/GO/GOMtMapped.txt", sep = '\t', stringsAsFactors = FALSE, header = FALSE)
colnames(GOF) <- c("Gene", "GO")
GOF <- cbind(GOF, value = 1)
xx <- as.list(GOBPANCESTOR)
xx <- xx[!is.na(xx)]
BIOP <- c(0)
for( i in 1:length(xx)){
  BIOP <- append(BIOP, names(xx[i]))}
#restrict GO terms that are only in a biological process
final_mat = acast(GOF, Gene ~ GO, fill=0, value.var="value")
#filter for genes only from the ones tested
test = final_mat[,colnames(final_mat) %in% BIOP,drop=F]
#command to run the enrichments
Results <- perform_enrichment(GeneMAT, test, minSize = 3, p.value.cutoff=0.05)
#add function column on the end
Results$Function=Term(as.vector(Results$id))
write.table(Results, "DEvsWTorNPD24.enrichment.txt", quote = FALSE, sep='\t')














#Organize --------------------
#AllDE = rownames(norm)
#AllDE = scan("~/Data/PDPs/DESeq/DElist.txt", what="", sep="\n")
#DEgenes = Sub.Ad2
#DEtable = norm[row.names(norm) %in% AllDE, ]
#ordercols <- colnames(DEtable)[c(1,2,3,13,14,15,4,5,6,16,17,18,7,8,9,19,20,21,10,11,12,22,23,24)]
#DEtable = DEtable[,match(ordercols, colnames(DEtable))]
ScaleDE = t(scale(t(norm)))
ScaleDE[is.na(ScaleDE)] <- 0
cor2 <- cor(ScaleDE, use="pairwise.complete.obs", method="pearson")
heatmap.2(cor2, na.rm=TRUE, dendrogram = 'none', Rowv = NULL, Colv=NULL, margin=c(8,8), key=FALSE,  trace='none')
#Pval < 0.05
Sub.Pv05.1 = subset(testNPD24, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2))$id
Sub.Pv05.2 = subset(testNPD245, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2))$id
Sub.Pv05.3 = subset(testNPD124, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2))$id
Sub.Pv05.4 = subset(testNPD24.NPD245, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2))$id
Sub.Pv05.5 = subset(testNPD24.NPD124, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2))$id
#write(Sub.Pv05.1, file = "Gm12DE.txt", ncolumns = 1, append = FALSE, sep = "\n")
#write(Sub.Pv05.2, file = "Gm85DE.txt", ncolumns = 1, append = FALSE, sep = "\n")
#write(Sub.Pv05.3, file = "Gm26DE.txt", ncolumns = 1, append = FALSE, sep = "\n")