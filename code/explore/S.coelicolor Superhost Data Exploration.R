
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
library(readr) #Need this library to load/read .csv files

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

strain <- factor(c(rep("WT",27),rep("MUT",27)),levels = c("WT","MUT")) 
x$samples$strain <- strain #should I create another DGE list with same data 
                           #and make the group be according to strain type instead of creating a new field?

#FILTERING STEP 

cpm <- cpm(x) 
lcpm <- cpm(x, log = T) #use logcpm for filtering, not normalization
keep.exprs <- rowSums(cpm > 1) >= 3 #Here we apply a filter that states that for each
                                    # gene at least 3 of the samples should have a CPM value higher than 1.

x <- x[keep.exprs,, keep.lib.sizes = FALSE] #Consider changing this filtering step to account for fact that
dim(x)                                      #we know that gene expression in the MUT will be low for the
                                            #deleted genes. ie: filter based on WT and non-deleted MUT genes

#NORMALIZATION USING TMM

x <- calcNormFactors(x, method = "TMM")
lcpmTMM = cpm(x,normalized.lib.sizes= TRUE, log=TRUE) #not necessary

# Visualize effect of TMM, so do I need to use CPM to incorporate TMM?? I thought it was an alternative
lcpm2TMM = cpm(x,normalized.lib.sizes= TRUE, log=TRUE)
lcpm2 = cpm(x,normalized.lib.sizes= FALSE, log=TRUE)
View(lcpm2)
View(lcpm2TMM)
boxplot(lcpm2)
title(main = "logCPM")
boxplot(lcpm2TMM)
title(main = "logCPM w/TMM")
remove(lcpm2)
remove(lcpm2TMM)

# PCA -ED

# Good to see if there are any outliers in the data: cant really tell from plot A tho
par(mfrow = c(1, 2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels = group, col = col.group)   #is multidimensional scaling = PCA??
title(main = "A. Color based on timegroup")
col.group = c(rep('#E41A1C',27),rep('#377EB8',27))
plotMDS(lcpm, labels = group, col = col.group, dim = c(1, 2)) 
title(main = "B. Color based on strain")

# Look interactively with Glimma
glMDSPlot(lcpm, labels = strain, groups = x$samples$strain,
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
s.class(pcar$li[,1:2], strain, cpoint = 1, col = c('blue','red'))

####################################################################
# Linear models

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
      
      # design <- model.matrix(~ strain * t1 * t2 * t3 * t4 * t5 * t6 * t7 * t8 * t9)

t1 <- rep(c(1,0,0,0,0,0,0,0,0),6)
t2 <- rep(c(0,1,0,0,0,0,0,0,0),6)
t3 <- rep(c(0,0,1,0,0,0,0,0,0),6)
t4 <- rep(c(0,0,0,1,0,0,0,0,0),6)
t5 <- rep(c(0,0,0,0,1,0,0,0,0),6)
t6 <- rep(c(0,0,0,0,0,1,0,0,0),6)
t7 <- rep(c(0,0,0,0,0,0,1,0,0),6)
t8 <- rep(c(0,0,0,0,0,0,0,1,0),6)
t9 <- rep(c(0,0,0,0,0,0,0,0,1),6)
design <- model.matrix(~ strain * ( t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9)) #i dont understand what all the terms in the design matrix mean
                                                          #when i include all of the time terms, specifically the ones with ":" in their name
design <- model.matrix(~ strain + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9) 

design <- model.matrix(~ strain ) #do i need to make a design model for each different time term or can i include them all in one?
design

x <- estimateDisp(x, design) #give me warning when i have all time interaction terms in design matrix

# Fit data to linear model
mfit <- glmQLFit(x, design)
# Inspect the coefficients that are found
head(mfit$coefficients)

# Perform DE analysis per factor
de <- glmLRT(mfit,coef = 2)
topTags(de)
tt <- topTags(de, n = Inf)
genes <- tt$table$genes[tt$table$FDR < 0.001]
plotMDS(lcpm[genes, ], labels = group, col = col.group) #just using the model with strain type, this plot looks exactly like the PCA plot 
                                    #previously obtained, so what was the point of going through the linear model?

#look at interaction terms

#get uniprot gs list
##################################################################

#benmodel

v <- voom(x,design,plot=T)#warning message? gives warning when more than one time term included in model
#also does VOOM do the same thing as estimateDisp?

fit <- lmFit(v,design) #this takes forever when i include all time terms in model
eb  <- eBayes(fit)
topTable(eb,coef = 2) #what does coef refer to? column of design matrix, we are interested in the ones after the intercepet right?