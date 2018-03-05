
## S.coelicolor Superhost Data Exploration
#  Author: Francisco Zorrilla

# PACKAGES
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "edgeR", "piano", "biomaRt", "Glimma"))
install.packages(c('tidyverse', 'RColorBrewer','snowfall','tibble'))

library(limma)
library(edgeR)
library(tidyverse) # Collection of useful R-packages
library(RColorBrewer)
library(biomaRt)
library(Glimma)
library(piano)
library(snowfall)
library(readr) #Need this library to load/read .csv files
library(tibble)

#CLEAR ENVIRONMENT

rm(list = ls())

#IMPORT DATA

setwd("C:/Users/zorrilla/Documents/GitHub/S.coelicolor/data/raw_internal/RNAseq") #Go to data folder

rawCountData_M1152 <- read_csv("rawCountData_M1152.csv") #These two lines give an error/warning message about a missing
rawCountData_M145 <- read_csv("rawCountData_M145.csv")   #column name, but they seem to work fine.

setwd("C:/Users/zorrilla/Documents/GitHub/S.coelicolor/code/explore") #Go back to code explore folder

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
remove(rawCountData_M1152,rawCountData_M145)

View(rawCounts) #Check that new data frame looks good

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
strain <- factor(c(rep("WT",27),rep("MUT",27)),levels = c("WT","MUT")) 
x$samples$strain <- strain 

#NORMALIZATION USING TMM

x <- calcNormFactors(x, method = "TMM") #do before filtering, violates assumption that downreg genes=upreg genes
#lcpmTMM = cpm(x,normalized.lib.sizes= TRUE, log=TRUE) #not necessary, non standard, need to justify reason for donig this


#FILTERING STEP 

cpm <- cpm(x,normalized.lib.sizes= TRUE) #i think cpm() uses norm.lib.sizes by default
#cpm2 <- cpm(x) #should there be a diffrence between cpm and cpm2?

#lcpm <- cpm(x,normalized.lib.sizes= TRUE, log=TRUE) #use logcpm for filtering, not normalization. Should i use this?

keep.exprs <- rowSums(cpm > 1) >= 3 #Here we apply a filter that states that for each
                                    # gene at least 3 of the samples should have a CPM value higher than 1.

        #is it better to filter using lcpm or cpm?

x <- x[keep.exprs,, keep.lib.sizes = FALSE] #Consider changing this filtering step to account for fact that
dim(x)                                      #we know that gene expression in the MUT will be low for the
                                            #deleted genes. ie: filter based on WT and non-deleted MUT genes


# Visualize effect of TMM, so do I need to use CPM to incorporate TMM?? I thought it was an alternative
lcpmTMM = cpm(x,normalized.lib.sizes= TRUE, log=TRUE)
lcpm2 = cpm(x,normalized.lib.sizes= FALSE, log=TRUE)
boxplot(lcpm2)
title(main = "logCPM")
boxplot(lcpmTMM)
title(main = "logCPM w/TMM")
remove(lcpm2)

#note, if using CPM or logCPM make sure it is recalculated after filtering

# PCA -ED

# Good to see if there are any outliers in the data: cant really tell from plot A tho
par(mfrow = c(1, 2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
x11();plotMDS(lcpmTMM, labels = group, col = col.group)   #use TMM not lcmp, use spearman
title(main = "A. Color based on timegroup")
col.group = c(rep('#E41A1C',27),rep('#377EB8',27))
x11();plotMDS(lcpmTMM, labels = group, col = col.group, dim = c(1, 2)) 
title(main = "B. Color based on strain")
            #one of the time point 9 samples from the mutant strain looks like it may be an outlier
            #especially considering the general "walk" pattern we see where the smaller numbers tend to be
            #lower than the bigger numbers

# Look interactively with Glimma
glMDSPlot(lcpmTMM, labels = strain, groups = x$samples$strain,
          launch = T) #plot is good but table is kind of wonky/doesnt display information besides strain

# PCA - Benj

require(ade4)

#Transpose the data
countsForPCA = t(lcpmTMM)

#Perform PCA on the counts
pcar <- dudi.pca(countsForPCA, center = TRUE, scale = FALSE, scannf = FALSE, nf=10)

#Check how much of the total variance each principal component accounts for:
var <- pcar$eig/sum(pcar$eig)*100
plot(var, type = 'b')

#Plot factorial map with representation of observations in the 1st 2 components:
x11();s.class(pcar$li[,1:2], strain, cpoint = 1, col = c('blue','red'))
s.class(pcar$li[,2:3], strain, cpoint = 1, col = c('blue','red'))


# LINEAR MODELS

      # strainTest <- c(rep(1,3),rep(0,3))
      # t1 <- rep(c(1,0,0,0,0,0,0,0,0),6)
      # t2 <- rep(c(0,1,0,0,0,0,0,0,0),6)
      # t3 <- rep(c(0,0,1,0,0,0,0,0,0),6)
      # t4 <- rep(c(0,0,0,1,0,0,0,0,0),6)
      # t5 <- rep(c(0,0,0,0,1,0,0,0,0),6)
      # t6 <- rep(c(0,0,0,0,0,1,0,0,0),6)
      # t7 <- rep(c(0,0,0,0,0,0,1,0,0),6)
      # t8 <- rep(c(0,0,0,0,0,0,0,1,0),6)
      # t9 <- rep(c(0,0,0,0,0,0,0,0,1),6)
      
      # design <- model.matrix(~ strain ) #only strain, no time component 
      # design <- model.matrix(~ strain + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9) #no interaction
      # design <- model.matrix(~ strain * t1 * t2 * t3 * t4 * t5 * t6 * t7 * t8 * t9) #dont use this, it creates
                                                                                      #weird interaction terms

t1 <- rep(c(1,0,0,0,0,0,0,0,0),6)
t2 <- rep(c(0,1,0,0,0,0,0,0,0),6)
t3 <- rep(c(0,0,1,0,0,0,0,0,0),6)
t4 <- rep(c(0,0,0,1,0,0,0,0,0),6)
t5 <- rep(c(0,0,0,0,1,0,0,0,0),6)
t6 <- rep(c(0,0,0,0,0,1,0,0,0),6)
t7 <- rep(c(0,0,0,0,0,0,1,0,0),6)
t8 <- rep(c(0,0,0,0,0,0,0,1,0),6)
t9 <- rep(c(0,0,0,0,0,0,0,0,1),6)
design <- model.matrix(~ strain * ( t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9)) 
design

#daniel: try making design matrix with one term for time, with differnt levels
  #if we do this then we cant get interaction terms right?

x <- estimateDisp(x, design) # Calculate dispersion using x and design model
mfit <- glmQLFit(x, design) # Fit data to linear model
head(mfit$coefficients) # Inspect the coefficients that are found

# Perform DE analysis per factor
de <- glmLRT(mfit,coef = 2) #is this legal to put groups of coefficients at the same time? or do i have to do 1 at a time
                                #how does choice of coeffcients influence pvalues?

topTags(de) #is the fold change relative to MUT or WT?
tt <- topTags(de, n = Inf) #store top DE genes in list tt
genes <- tt$table$genes[tt$table$FDR < 0.001]

#plotMDS(lcpmTMM[genes, ], labels = group, col = col.group) #not needed

#look at interaction terms

################################do tSNE

library(Rtsne)


##### GSA

## Download GO terms from Uniprot, S. coelicolor 'proteome': http://www.uniprot.org/uniprot/?query=proteome:UP000001973
## Use columns 'Gene names (ordered locus)' and 'Gene ontology (GO)'
## Save as tab-delimited
GO<-read.delim('uniprot-proteome%3AUP000001973.tab')
GO<-GO[,-1] # Remove 'Entry' column
GO<-GO[!GO$Gene.names...ordered.locus..=='',] # Remove rows without gene name
colnames(GO)<-c('gene','GO') # Rename columns

# install.packages('splitstackshape') # Run this installation command if the package is not installed yet.
library(splitstackshape)
GO<-cSplit(GO,"GO",sep=';',direction='long') # Split the GO column, separated by ';', and make a long list
head(GO)

library(piano)
GO<-loadGSC(GO) #ignore warning?

fc <- subset(tt$table,select = "logFC")
p <- subset(tt$table,select = "FDR")

##############BEN
resGSA  <- runGSA(geneLevelStats=p, #We supply our p-values
                  directions=fc,    #We supply the fold-changes
                  geneSetStat="reporter", #We decided the method
                  signifMethod="geneSampling", #We choose the statistical distribution against which we measure our confidence
                  adjMethod="fdr",  #We adjust the gene-set p-values for multiple testing
                  gsc=GO, #We supply which gene-sets to test
                  ncpus=4) #We decide the number of processor cores to run the calculations

windows()
GSAheatmap(resGSA,cutoff =10,adjusted = T,ncharLabel = Inf, cellnote='pvalue') #cutoff of 9 - 10 is readable
#dev.off() 

#################ED
gsaRes <- runGSA(p, fc, gsc = GO, gsSizeLim = c(10, 400))

x11();GSAheatmap(gsaRes, cutoff = 5, cellnote = "nGenes",ncharLabel = Inf)

#BEN/ED GSAs give me different heatmaps: numbers are different(FCs?) and the GO terms
#are also somewhat differrent. which one should i use and what parameters are causing
#these differences? Why did Ed narrow down the gsSizeLim to c(10,400)?
#How should I determine the parameters to set for geneSetStat and signifMethod?

#how should i run GSA in order to determine DE expression between time points?
#start from line de <- glmLRT(mfit,coef = 2) but change coef to an interaction term coef