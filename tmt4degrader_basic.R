#######
# Owner: Fischer lab
# Author: Eric S. Fischer
# Last modified: ESF 2018-09-23
# Dependencies: Fischer lab helper functions
# R including packages: limma, ggplot
#######

# Define input files and variables
#rawdata file: This should be in typical format providing RAW reporter ion intensities (not scaled/normalized).
#columns should be: Protein ID; Gene Symbol; Description; Quantified peptides; TMT channels 1-10 reporter ion intensities
#of advantage to provide meaningful names for TMT channels in input file
infile <- "rawdata/experiment.txt"
#number of minimal spectral counts required to use data
sc <- 2
#dataset prefix for naming of output files
pre <- "esf-wp001"
# use intensities (1) or spectral counts (2) as weights for limma
w <- 1
# Define compound to look at - based on colnames - Contrast is formed between the chosen compound, and the others are used as controls.
# DMSO is identified based on colname and used as control
contrast <- "ESF-01"
# Number of TMT channels used:
tmt <- 10


# load dependencies
library(limma)



### load TMT data - NEEDS NO STANDARD INPUT FILES MAY NEED ADJUSTMENT TO FIT DATA CONVENTIONS ####
D.exp <- read.delim(infile, header=TRUE, sep="\t", as.is=TRUE)
dim(D.exp) # [1] 9931 14
colnames(D.exp)
D <- D.exp
colnames(D)
summary(D)
# Make data uniform
# Protein IDs, Gene Symbol, Description, Quantified Spectral Counts, TMT channels
D <- D[,c(1,2,3,4,5:(5+tmt-1))]
colnames(D)
colnames(D)[1:4] <- c("ProteinID.Uniprot","Gene.Symbol","Description","Quantified.peptides")
colnames(D)

####################### STANDARD SCRIPT FROM HERE ON ####################
# remove data with spectral count < min spectral counts defined above
sel1 <- D$Quantified.peptides >= sc #& rowSums(D[,c(9:10)]) >= 10^2
summary(sel1)
D <- subset(D, sel1)

# normalize data
x <- colSums(D[,c(5:(5+tmt-1))]) #calculate sum of intensities for each channel
x <- x/max(x) # calculate normalization factors
x
y <- t(apply(D[,c(5:(5+tmt-1))],x, "/", x)) # apply normalization factos
Dn <- D
Dn[,c(5:(5+tmt-1))] <- y # Dn contains normalized value
rm(x); rm(y)
summary(Dn)

# create matrix for abundance
A <- as.matrix(Dn[,c(5:(5+tmt-1))])
colnames(A)
rownames(A) <- Dn$Gene.Symbol
summary(A)

# Scaling protein TMT channels to 100
x <- rowSums(A)
ZZ <- A/x
zz <- ZZ*100
# create matrix for relative abundance
a <- zz
rm(zz); rm(ZZ); rm(x)
a <- log2(a)

### Define experimental factors
treat <- factor(ifelse(grepl(contrast, colnames(a)),"DRUG","CTRL"), levels=c("CTRL","DRUG"))
ctrl <- factor(ifelse(grepl("DMSO",colnames(a)),"DMSO","TREAT"), levels=c("DMSO", "TREAT"))
treat
ctrl


######## LIMMA MODERATED T-STATISTICS #####
design <- model.matrix(~ treat + ctrl) # define matrix
rownames(design) <- colnames(a) # assign names

design
### matrix must be of full rank: Matrix::rankMatrix(design)==ncol(design)
# rund ebayes fit with limma
if(w == 1){
  fit <- lmFit(a, design=design, weights=rowSums(A)) # with intensities as weight
} else {
  fit <- lmFit(a, design=design, weights=D$Quantified.spectral.counts) # with spectral counts as weight
}
fit <- eBayes(fit)
tt <- topTable(fit, coef="treatDRUG", genelist=rownames(a), number=nrow(a))
options(width=160); head(tt, n=50)
#Write hit table
write.csv(tt, file=paste("tables/",pre, "_Limma_output_DRUG_vs_CTRLs.csv"))


# Plot & Hits by adjusted p.value
png(filename = paste("plots/",pre,"_hits_adjPvalue.png"), width = 1200, height = 1200)
sel <- tt$logFC <= -log2(1.25) & tt$adj.P.Val <= 0.2 | tt$logFC >= log2(1.25) & tt$adj.P.Val <= 0.2
summary(sel)
sel2 <- tt$logFC <= -log2(1.25) & tt$adj.P.Val <= 0.5 & tt$adj.P.Val >= 0.05
summary(sel2)
#plot
par(mar=c(5,4,4-3,2)+.1)
plot(log10(tt$adj.P.Val),tt$logFC,  xlab="-log10 adjusted P value", ylab="logFC (log2 TREAT / CTRL)",
     #ylim=c(-2.2,2.2), 
     #xlim=c(-2.2,2.2),
     cex=0.8,
     pch=20, col=ifelse(sel, "red", "#22222222"))
text(log10(tt[sel,"adj.P.Val"]), tt[sel,"logFC"],adj = c(-0.25,0.55), labels=tt[sel,"ID"], col="black", cex=1.0)
abline(v=log10(0.05), lty=2)       # P value < 0.05
abline(h=c(-1,1)*log2(1.25), lty=2) # > 25% change
legend(pch=20, "topleft", bty="n", legend="Hits: adj p < 0.05; logFC > 25%", cex=1.0, col="red")
text(-2.5, 0.4, labels="25% up-regulated")
text(-2.4, -0.45, labels="25% down-regulated")
text(-1,-2.7, labels = "p = 0.05")
dev.off()


# Plot & Hits by p.value
pdf(file=paste("plots/",pre,"_hits_Pvalue.pdf"), width = 16, height = 16, useDingbats = FALSE)
sel <- tt$logFC <= -log2(1.25) & tt$P.Value <= 0.01 | tt$logFC >= log2(1.25) & tt$P.Value <= 0.01
summary(sel)
#plot
par(mar=c(5,4,4-3,2)+.1)
plot(log10(tt$P.Value),tt$logFC,  xlab="-log10 P value", ylab="logFC (log2 TREAT / CTRL)",
     ylim=c(-1.3,1.3), 
     cex=0.8,
     pch=20, col=ifelse(sel, "red", "#22222222"))
text(log10(tt[sel,"P.Value"]), tt[sel,"logFC"],adj = c(-0.25,0.55), labels=tt[sel,"ID"], col="black", cex=1.0)
abline(v=log10(0.01), lty=2)       # P value < 0.01
abline(h=c(-1,1)*log2(1.25), lty=2) # > 50% change
legend(pch=20, "topleft", bty="n", legend="Hits: p < 0.01; logFC > 25%", cex=1.0, col="red")
text(-6.5, 0.41, labels="25% up-regulated")
text(-6.4, -0.45, labels="25% down-regulated")
text(-1.4,-1.3, labels = "p = 0.01")
dev.off()