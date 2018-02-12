
## S.coelicolor Superhost Data Exploration
#  Author: Francisco Zorrilla

# PACKAGES
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "edgeR", "piano", "biomaRt", "Glimma"))
install.packages(c('tidyverse', 'RColorBrewer'))

# After installing the software the first time, load the R-packages that we need. This you need
# to do every time you start a new R session.
library(limma)
library(edgeR)
library(tidyverse) # Collection of useful R-packages
library(RColorBrewer)
library(biomaRt)
library(Glimma)
library(piano)
library(snowfall)

#CLEAR ENVIRONMENT

rm(list = ls())

#IMPORT DATA

setwd("~/Superhost Individual Project/data/raw_internal/RNAseq") #Go to data folder

library(readr) #Need this library to load/read .csv files
rawCountData_M1152 <- read_csv("rawCountData_M1152.csv") #These two lines give an error/warning message about a missing
rawCountData_M145 <- read_csv("rawCountData_M145.csv")   #column name, but they seem to work fine.

setwd("~/Superhost Individual Project/code/explore") #Go back to code explore folder

#CHECK DATA

head(rawCountData_M1152) #Just to make sure data was properly imported/nothing weird happened
head(rawCountData_M145)

dim(rawCountData_M1152) #Check dimensions of count matrices, should be equal
dim(rawCountData_M145)

#MAKE SOME COSMETIC MODIFICATIONS

row.names(rawCountData_M1152)=rawCountData_M1152$X1 #Set row names to gene names. These two lines give a warning
row.names(rawCountData_M145)=rawCountData_M145$X1   #message about tibble depreciation?? They work fine.

rawCountData_M1152$X1 <- NULL # delete columns with gene names
rawCountData_M145$X1 <- NULL

rawCounts = cbind(rawCountData_M145,rawCountData_M1152) #Create single dataframe combining raw counts
                                                        #of the two strains

view(rawCounts) #Check that new data frame looks good

    #Columns 1-9 contain 9 time points P(6,14,18,22,26,30,34,38,42) of F516 sample (WT)
    #Columns 10-18 contain 9 time points P(6,14,18,22,26,30,34,38,42) of F517 sample (WT)
    #Columns 19-27 contain 9 time points P(6,14,18,22,26,30,34,38,42) of F518 sample (WT)
    #Columns 28-36 contain 9 time points P(18,26,30,34,38,42,46,50,51) of F519 sample (MUT)
    #Columns 37-45 contain 9 time points P(18,26,30,34,38,42,46,50,51) of F521 sample (MUT)
    #Columns 46-54 contain 9 time points P(18,26,30,34,38,42,46,50,51) of F522 sample (MUT)

# Move data into DGEList, which can store additional metadata
x <- DGEList(counts = rawCounts, genes = rownames(rawCounts))


# Add the grouping information (and enforce the order of the samples).


group <- factor(rep(1:9,6),levels = c("1","2","3","4","5","6","7","8","9"))
x$samples$group <- group

#FILTERING STEP 

cpm <- cpm(x) 
lcpm <- cpm(x, log = T) #use logcpm for filtering, not normalization
keep.exprs <- rowSums(cpm > 1) >= 3 #Here we apply a filter that states that for each
                                    # gene at least 3 of the samples should have a CPM value higher than 1.
x <- x[keep.exprs,, keep.lib.sizes = FALSE]
dim(x)

#NORMALIZATION USING TMM

# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x <- calcNormFactors(x, method = "TMM")
# If we look at the normalization factors, we see that the effect of TMM normalization was
# quite mild, which might not be too surprising with similar library sizes for each sample
x$samples

#================== 4. Unsupervised clustering =======EDIT THIS SECTION====================

# Good to see if there are any outliers in the data:
par(mfrow = c(1, 2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1") #set colors to WT, look at col.group
col.group <- as.character(col.group)
plotMDS(lcpm, labels = group, col = col.group)
title(main = "A. Dimensions 1&2")
plotMDS(lcpm, labels = group, col = col.group, dim = c(3, 4))
title(main = "B. Dimensions 3&4")

# This gives question: is one of the ref and one of the pH swapped?
# Something that can explain separation in dim 3? (same reactors?)

# Look interactively with Glimma
glMDSPlot(lcpm, labels = group, groups = x$samples[, 2], folder = "scratch/glimma",
          launch = T)

# Let's leave all samples in for now, and look what the rest of the analysis has to offer.

#DATA EXPLORATION, ADD BOXPLOTS COMPARING TIME POINTS

boxplot(logCPM[,],ylab="logCPM") #Should this be centered around zero like in the bioinformatics ex.3?

boxplot(log(filtCounts+1),ylab="log(rawgeneExp+1)") #Comparing this to the logCPM plot shows effect of normalization

rawChangeInExpression =rowSums(rawCountData_M145) - rowSums(rawCountData_M1152) #look at min and max of raw
which(rawChangeInExpression==max(rawChangeInExpression))

#interesting gene (SCO6282) Identified by rawChangeInExpression plot above CHANGED AFTER FILTRATION
intGene1WT = t(logCPM["SCO6282",1:27])
intGene1MUT = t(logCPM["SCO6282",28:54])
intGene1logCPM = cbind(intGene1WT,intGene1MUT)
boxplot(intGene1logCPM, ylab = 'logCPM(Gene Expression)') #still differentially expressed after normalization
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3293853/ ----PAPER THAT MENTIONS THIS GENE/ENZYME

which(rawChangeInExpression==min(rawChangeInExpression)) #find most upregulated in MUT
intGene2WT = t(logCPM["SCO4296",1:27])
intGene2MUT = t(logCPM["SCO4296",28:54])
intGene2logCPM = cbind(intGene2WT,intGene2MUT)
boxplot(intGene2logCPM, ylab = 'logCPM(Gene Expression)') #not very differentially expressed after normalization

#LINEAR MODEL TO DETERMINE DE OF GENE SCO6282 & SCO4296

type = c(rep("WT",27), rep("MUT",27)) #create type variable for distinguishing between WT and MUT
type = as.factor(type)
type = relevel(type, ref = "WT")

SCO6282data=cbind(type,t(logCPM["SCO6282",]))
SCO6282data=data.frame(SCO6282data)
fit1= lm(SCO6282~type, data=SCO6282data)

SCO4296data=cbind(type,t(logCPM["SCO4296",]))
SCO4296data=data.frame(SCO4296data)
fit2= lm(SCO4296~type, data=SCO4296data)

#LINEAR MODEL FOR ALL GENES

pvalsFit1 = c()
estimatesFit1 = c()

for (z in 1:7852){ #use edgeR and limma
  geneData = cbind(type,t(logCPM[z,]))
  geneData = data.frame(geneData)
  fit1X = lm(geneData[,2]~type, data=geneData)
  estimatesFit1[z] = summary(fit1X)$coeff[2]
  pvalsFit1[z] = summary(fit1X)$coeff[8]
}

sigGeneNames = row.names(logCPM)
fitDataFrame = cbind(pvalsFit1,estimatesFit1)
fitDataFrame = data.frame(fitDataFrame)
row.names(fitDataFrame) = sigGeneNames

pvalsFit1Adj = p.adjust(fitDataFrame[,1], method = "fdr", n = 7852)
fitDataFrame = cbind(fitDataFrame,pvalsFit1Adj)

sigPvalsFit1idx = pvalsFit1Adj<0.05
#sigPvalsFit1 = pvalsFit1Adj[sigPvalsFit1idx] CAN DELETE
sigPvalsFit1DataFrame = fitDataFrame[fitDataFrame[,3]<0.05,]

#sigEsimates1 = fitDataFrame[sigPvalsFit1idx,2] CAN DELETE

#sigGeneFit1DataFrame = cbind(sigPvalsFit1,sigEsimates1) CAN DELETE
#sigGeneFit1DataFrame = data.frame(sigGeneFit1DataFrame) CAN DELETE

fit1up = sigPvalsFit1DataFrame[sigPvalsFit1DataFrame[,2]>0,]
fit1down = sigPvalsFit1DataFrame[sigPvalsFit1DataFrame[,2]<0,]

#HEAT MAP OF TOP 100 SIG DE GENES

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

sortedPvals = sort(fitDataFrame[,3],decreasing=FALSE)
hunnedCutoffIdx = sortedPvals[100]
hunnedFrame = fitDataFrame[fitDataFrame[,3]<=hunnedCutoffIdx,]
hunnedGenes = row.names(hunnedFrame)
tlogCPM = t(logCPM)
xsig = tlogCPM[,hunnedGenes]
library(RColorBrewer)
library(gplots)
mypalette = brewer.pal(11,"RdYlBu")
morecols = colorRampPalette(mypalette)
mycols=rev(morecols(255))
column.cols=c("purple","orange")[type]
pdf("top100sigGenesHeatmap.pdf",height=9,width=6)
heatmap.2(t(xsig),trace='none',col=mycols,main='The 100 most significant
          genes',ColSideColors=column.cols)
dev.off()

#LINEAR MODEL BY TIME POINT

#Columns 1-9 contain 9 time points P(6,14,18,22,26,30,34,38,42) of F516 sample (WT)
#Columns 10-18 contain 9 time points P(6,14,18,22,26,30,34,38,42) of F517 sample (WT)
#Columns 19-27 contain 9 time points P(6,14,18,22,26,30,34,38,42) of F518 sample (WT)
#Columns 28-36 contain 9 time points P(18,26,30,34,38,42,46,50,51) of F519 sample (MUT)
#Columns 37-45 contain 9 time points P(18,26,30,34,38,42,46,50,51) of F521 sample (MUT)
#Columns 46-54 contain 9 time points P(18,26,30,34,38,42,46,50,51) of F522 sample (MUT)

timePt = c(rep("1",3),rep("2",3),rep("3",3),rep("4",3),rep("5",3),rep("6",3),
           rep("7",3),rep("8",3),rep("9",3),rep("1",3),rep("2",3),rep("3",3),
           rep("4",3),rep("5",3),rep("6",3),rep("7",3),rep("8",3),rep("9",3))


#timePt = c(rep("P6",3),rep("P14",3),rep("P18",3),rep("P22",3),rep("P26",3),rep("P30",3),
 #          rep("P34",3),rep("P38",3),rep("P42",3),rep("P18",3),rep("P26",3),rep("P30",3),
  #         rep("P34",3),rep("P38",3),rep("P42",3),rep("P46",3),rep("P50",3),rep("P51",3))



typePvalsFit2 = c()
timePtPvalsFit2 = c()
typeEstimatesFit2 = c()
timePtEstimatesFit2 = c()
timePt =as.factor(timePt)

for (z in 1:7852){
  geneData = cbind(type,timePt, t(logCPM[z,]))
  geneData = data.frame(geneData)
  fit2X = lm(geneData[,3]~type+timePt, data=geneData)
  typePvalsFit2 = c()
  timePtPvalsFit2 = c()
  typeEstimatesFit2 = c()
  timePtEstimatesFit2 = c()
  estimatesFit1[z] = summary(fit1X)$coeff[2]
  pvalsFit1[z] = summary(fit1X)$coeff[8]
}

#PCA

pca=prcomp(tlogCPM)
summary(pca)

plot(pca$x[,1],pca$x[,2],col =column.cols, xlab = "PC1",ylab = "PC2")
title("PCA Scatter Plot")

plot(pca$x[,1],pca$x[,3],col =column.cols, xlab = "PC1",ylab = "PC3")
title("PCA Scatter Plot")

plot(pca$x[,2],pca$x[,3],col =column.cols, xlab = "PC2",ylab = "PC3")
title("PCA Scatter Plot")
