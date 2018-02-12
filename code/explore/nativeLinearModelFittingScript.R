#OLD CODE FOR LINEAR MODEL FITTING + PCA, GONNA USE PACKAGES INSTEAD OF USING SUCH A NATIVE APPROACH,
#SHOULD STILL BE A GOOD REFERENCE. COMPARE WITH RESULTS OBTAINED FROM PACKAGES LATER


#LINEAR MODEL TO DETERMINE DE OF GENE SCO6282 & SCO4296

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
