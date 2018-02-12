# R-script accompanying CHASSY transcriptomics workshop BASED SCRIPT

# Load R-packages required for analysis
# First time: need to install the packages:

# In case you don't have admin rights on your computer, make sure that
# R-packages are installed in a folder of choice. Adjust and uncomment
# the next line and run every time you start a new R-session:
# .libPaths(c("C:/Box Sync/Documents/Work/Teaching/Chassy_omics/scratch/Rlib",.libPaths()))

# There are two major sources for R packages, CRAN (https://cran.r-project.org/) and
# Bioconductor (https://bioconductor.org/). CRAN packages tend to be more general, while
# Bioconductor is biology-related. Installation commands are slightly different.
# These packages are from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "edgeR", "piano", "biomaRt", "Glimma"))
# And these are from CRAN
install.packages(c('tidyverse', 'RColorBrewer'))

# Once installed, you don't need to repeat these commands above (you can comment them out with
# a #-sign in front).

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


#================== 1. Loading and restructure data ================================================

# Import the raw counts as provided by featureCounts.
# Skip first line, contains parameters from featureCounts
# counts <- read.delim('rawCountData_M145.csv', skip = 1)

library(readr)
rawCountData_M1152 <- read_csv("rawCountData_M1152.csv") #These two lines give an error/warning message about a missing
rawCountData_M145 <- read_csv("rawCountData_M145.csv")   #column name, but they seem to work fine.
row.names(rawCountData_M1152)=rawCountData_M1152$X1 #Set row names to gene names. These two lines give a warning
row.names(rawCountData_M145)=rawCountData_M145$X1   #message about tibble depreciation?? They work fine.
rawCountData_M1152$X1 <- NULL # delete columns with gene names
rawCountData_M145$X1 <- NULL
rawCounts = cbind(rawCountData_M145,rawCountData_M1152) #Create single dataframe combining raw counts
#of the two strains

# Move data into DGEList, which can store additional metadata
x <- DGEList(counts = rawCounts, genes = rownames(rawCounts))


# Add the grouping information (and enforce the order of the samples).

group <- factor(c(rep("1",3),rep("2",3),rep("3",3),rep("4",3),rep("5",3),rep("6",3),
           rep("7",3),rep("8",3),rep("9",3),rep("1",3),rep("2",3),rep("3",3),
           rep("4",3),rep("5",3),rep("6",3),rep("7",3),rep("8",3),rep("9",3)),
           levels = c("1","2","3","4","5","6","7","8","9"))


x$samples$group <- group

# Now inspect what the data looks like:
x

# The content can also be summarized with:
dim(x)
# Showing the number of genes and the number of samples

#================== 2 Filtering low reads ================================================

# In x$samples$lib.size we already see that the number of reads varies a bit.
# Let's normalize for library size, by taking the log count-per-million. (In lecture is
# discussed how (log)CPM is not a satisfactory normalization, but here we just apply it
# to filter low reads!).
cpm <- cpm(x)
lcpm <- cpm(x, log = T)

# With the following code we can make a plot showing the raw logCPM reads, the dotted
# line indicates zero logCPM (= 1 cpm)
nsamples <- ncol(x)
col <-
  brewer.pal(nsamples, "Paired") # Ignore warning, some samples will have identical color.
par(mfrow = c(1, 2))
plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.35),
  las = 2,
  main = "A. Raw data",
  xlab = "Log-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}

# There are many genes with low reads. Different filters can be used to
# remove these low counts. Here we apply a filter that states that for each
# gene at least 3 of the samples should have a CPM value higher than 1. For
# the smallest library (ana.5_S19) that means 11 reads.

keep.exprs <- rowSums(cpm > 1) >= 3
x <- x[keep.exprs,, keep.lib.sizes = FALSE]
dim(x)
# Now we have less genes (compared to when we ran dim(x) above)

# Alternatively, only keep those genes for which we have data in all samples:
x <- DGEList(counts = rawCounts, genes = rownames(rawCounts)) # Again, move raw counts into x
x$samples$group <- group
keep.exprs <- rowSums(cpm > 1) >= 15
x <- x[keep.exprs,, keep.lib.sizes = FALSE]
dim(x)

lcpm <- cpm(x, log = T)

plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.35),
  las = 2,
  main = "B. Filtered data",
  xlab = "Log-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}

#================== 3. Normalize counts ================================================

# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x <- calcNormFactors(x, method = "TMM")
# If we look at the normalization factors, we see that the effect of TMM normalization was
# quite mild, which might not be too surprising with similar library sizes for each sample
x$samples

# To visualize the effect of TMM, we make a fake data-set, where the counts from the first sample
# are reduced 20-fold, while the counts from the second sample are increased 5-fold.
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

# We then visualize this with two graphs, before and after normalization:
par(mfrow = c(1, 2))
lcpm2 <- cpm(x2, log = TRUE)
boxplot(lcpm2,  las = 2, col = col, main = "")
title(main = "A. Example: Unnormalised", ylab = "Log-cpm")
x2 <- calcNormFactors(x2)
lcpm2 <- cpm(x2, log = TRUE)
boxplot(lcpm2, las = 2, col = col, main = "")
title(main = "B. Example: Normalised", ylab = "Log-cpm")

# To avoid confusion, we'll remove x2 and lcmp2, as we only used this to demonstrate the TMM
# normalization in an extreme case
remove(x2)
remove(lcpm2)

#================== 4. Unsupervised clustering ================================================

# Good to see if there are any outliers in the data:
par(mfrow = c(1, 2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels = group, col = col.group)
title(main = "A. Dimensions 1&2")
plotMDS(lcpm, labels = group, col = col.group, dim = c(3, 4))
title(main = "B. Dimensions 3&4")

# This gives question: is one of the ref and one of the pH swapped?
# Something that can explain separation in dim 3? (same reactors?)

# Look interactively with Glimma
glMDSPlot(lcpm, labels = group, groups = x$samples[, 2],
          launch = T)

# Let's leave all samples in for now, and look what the rest of the analysis has to offer.

#======================= 5. Pairwise DE analysis ==============================================
# We can analyze using pairwise contrasts, or using a generalized linear model. Here is
# demonstrated how to do pairwise.
# Using the grouping, we can now do an additional correction for the mean-variance correlation 
x <- estimateDisp(x)
# We can plot the mean-variance correlation
par(mfrow = c(1, 1))
plotBCV(x)
# The tagwise (individual genes) displacement from the common variance will be taken into
# account in remaining statistical framework.

# Regardless, we can now get lists of differentially expressed genes for each comparison, using
# the labels specified in group (let's have a look again what we specified:)
group
de <- exactTest(x, pair = c('ref', 'anox')) # Reference first!
# This will compare 'ref' with 'anox'. Show top 20 results:
topTags(de, n = 20)

# Construct a table with all genes (n = Inf)
tt <- topTags(de, n = Inf)

# Write CSV file if you prefer to inspect in e.g. Excel. Note that we're
# writing it to the results folder.
write.csv(tt, file = "results/ref_anox_DEgenes.csv", row.names = F)

# Make a smear plot, highlighting as DE those genes that have FDR < 0.01
plotSmear(x, pair = c('ref', 'anox'), de.tags = rownames(
  tt$table)[tt$table$FDR < 0.01 & abs(tt$table$logFC) > 1])

# If you want to save a plot, you can easily do this by (1) specifying the pdf
# function, (2) making the plot, and (3) switching off the writing of files.
pdf("results/smearPlot_anoxVSref.pdf")
plotSmear(x, pair = c('ref', 'anox'), de.tags = rownames(
  tt$table)[tt$table$FDR < 0.01 & abs(tt$table$logFC) > 1])
dev.off()

# Make histogram of p-values, add vertical line at p = 0.05
hist(tt$table$PValue, breaks = 100)
abline(v = 0.05, col = 'red')

# Make volcano plot
volcanoData <- cbind(tt$table$logFC, -log10(tt$table$FDR))
colnames(volcanoData) <- c("logFC", "-log10Pval")
plot(volcanoData)
abline(v = c(-2, 2), col = 'red')
abline(h = 3, col = 'red')

# Adjust the code above to make tables and CSV files for all comparisons that you are
# interested in.

# Also an option to make interactive graph with glimma:
de <- exactTest(x) # Reference first!
dt <- decideTests(de)
glMDPlot(de, status = dt, counts = x$counts, groups = x$samples$group, transform = TRUE,
         folder = "scratch/glimma")

#================== 6. DE analysis using linear model===========================================
# Specify design matrix as detailed in lecture
design <- model.matrix(~group)
colnames(design) <- gsub("group", '', colnames(design))
design

# Again calculate dispersion, now further informed by design matrix
x <- estimateDisp(x, design)

# Fit data to linear model
mfit <- glmQLFit(x, design)
# Inspect the coefficients that are found
head(mfit$coefficients)

# Perform DE analysis per factor
de <- glmLRT(mfit, coef = 2)
topTags(de)

# Or all de <- glmLRT(fit, coef=c(2,3,4,5))
de <- glmLRT(mfit, coef = c(2,3,4,5))
topTags(de)
tt <- topTags(de, n = Inf)
write.csv(tt, file = "results/glm_DEgenes.csv", row.names = F)

# Make clustering plot of DE genes only (similar as shown for proteomics data)
# First get a list of genes that are strongly differentially expressed upon the conditions
genes <- tt$table$genes[tt$table$FDR < 0.001]
plotMDS(lcpm[genes, ], labels = group, col = col.group)

#================== 7. Gene-set analysis =========================================
# Ensembl (ensembl.org) is a database with a lot of genomics data. We will
# extract GO term annotation from S. cerevisiae from this database, using
# the biomaRt tool.
# Select ensembl database and S. cerevisiae dataset:
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "scerevisiae_gene_ensembl" ,host = "www.ensembl.org")
# Extract S. cerevisiae specific GO terms
mapGO <- getBM(attributes = c('ensembl_gene_id','go_id','name_1006','namespace_1003'), mart = ensembl)
# Remove blanks ("")
mapGO <- mapGO[mapGO[,2] != "",]
# Check the 10 first rows to see what we got:
mapGO[1:10,]

# Only look at biological process, first indicate what rows to keep
keepRows <- mapGO$namespace_1003 == 'biological_process'
# Then subset mapGO to only keep those rows
mapGO <- mapGO[keepRows, ]

# Only keep columns with gene and GO term description (1st and 3rd column)
mapGO <- mapGO[, c(1,3)]
head(mapGO)

# Instead of using functions such as bioMart to import GO terms, you can also
# make your own list of gene-sets and load them as CSV file
mapGO <- read.csv('data/ensemblGOterms.csv')

# Load gene-sets for usage by Piano
myGsc <- loadGSC(mapGO)

# GSA is always on one statistic, so for instance one comparison between
# 2 conditions. Here we extract one of the coefficients of the GLM, and its
# associated P-value. As detailed in the lecture, the coefficients of the GLM
# represents the log2FC caused the the associated factor.

# Here we look at the second factor, which is anoxia (see columnnames)
head(mfit$coefficients)

de <- glmLRT(mfit, coef = 2) # Fit linear model to extract the effect of anoxia
tt <- topTags(de, n = Inf) # Make a table with P-values and log2FC that will
# be used by Piano

# Extract necessary data from table
Pval <- tt$table$PValue
names(Pval) <- tt$table$genes
FC <- tt$table$logFC
names(FC) <- tt$table$genes

# Run GSA with the extracted P-value, FC, the previously loaded gene-set
# collection myGsc, and we set a size limit of 10, 400.
gsaRes <- runGSA(Pval, FC, gsc = myGsc, gsSizeLim = c(10, 400))

# We can show the results in a heatmap
GSAheatmap(gsaRes, cutoff = 1, cellnote = "nGenes")

# Sometimes R comes with an error message regarding plotting,
# for instance "Error in plot.new() : figure margins too large". If the
# picture is plotted as intended, you can just ignore these messages.
# If the picture is not shown, it can help to run the following command
# (potentially multiple times, until it returns "null device 1":

# dev.off()

# NADH oxidation might be interesting, let's have a closer look.
# We can find out which genes are part of this geneset:
# (make sure you use the correct quotation mark around `GO term`, this
# is needed due to spaces in the GO terms.)
genes <- gsaRes$gsc$`NADH oxidation`
# Extract the gene level data (note that this is for the anox coefficient)
tt[genes,]

# Maybe we want to see transcript levels of these genes for all samples.
# We'll use the normalized counts for this:
lcpm <- cpm(x)
boxData <- data.frame(lcpm[genes,], check.names = F)

# We make boxplots with ggplot2, which requires a specific data format (compare boxData
# before and after the following commands)
boxData <- cbind(genes = rownames(boxData), boxData) # Add genes as separate column
boxData <- gather(boxData, "sample", "count", 2:16) # Reorganize structure for ggplot2
map <- setNames(group, colnames(lcpm)) # Replicates should have same sample name. Map them here,...
boxData$sample <- map[boxData$sample] # and replace replicate name with sample name here.

# Make boxplots with ggplot2
ggplot(boxData, aes(sample, count)) +
  geom_boxplot(aes(colour = sample), size = 1) +
  facet_wrap( ~ genes, scales = "free_y") +  theme_bw() +
  labs(title = "NADH oxidation", x = "Sample", y = "Normalized counts")

# We can also connect the GO terms in a network plot
networkPlot(gsaRes, 'distinct', 'down',
            adjusted = T, significance = 0.01)

# Or show the results in a table
gsaTable <- GSAsummaryTable(gsaRes)
View(gsaTable)

# Run hypergeometric GSA, where we manually select a set of genes. Here, we
# select the most DE genes (FDR < 0.01) in high temperature (coef 5).
de <- glmLRT(mfit, coef = 5)
tt <- topTags(de, n = Inf)

hypGsaRes <- runGSAhyper(genes = tt$table$genes, pvalues = tt$table$FDR,
                         pcutoff = 0.01, gsc = myGsc)

networkPlot(hypGsaRes, 'non')

# Run consensus GSA. Let's speed up the process, we will only do 1000
# permutations when calculating the gene-set signficance (default is 10,000).
# To speed things up, we can check how many CPUs your computer has, and
# try to use them all
library(parallel)
cores <- as.numeric(detectCores())
# This is the nubmer of cores available:
cores
# Okay, back to consensus GSA:

gsaRes1 <- runGSA(Pval, directions = FC, geneSetStat = "mean", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes2 <- runGSA(Pval, directions = FC, geneSetStat = "median", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes3 <- runGSA(Pval, directions = FC, geneSetStat = "sum", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes4 <- runGSA(FC, geneSetStat = "maxmean", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes5 <- runGSA(Pval, directions = FC, geneSetStat = "fisher", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes6 <- runGSA(Pval, directions = FC, geneSetStat = "stouffer", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes7 <- runGSA(Pval, directions = FC, geneSetStat = "tailStrength", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes8 <- runGSA(FC, geneSetStat = "gsea", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes9 <- runGSA(FC, geneSetStat = "page", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)

# Once we have several GSA results, we can combine them and make a heatmap
resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5, gsaRes6,
                gsaRes7, gsaRes8, gsaRes9)
names(resList) <- c("mean", "median", "sum", "maxmean", "fisher", "stouffer",
                    "tailStrength", "gsea", "page")
ch <- consensusHeatmap(resList, cutoff = 5, method = "mean")