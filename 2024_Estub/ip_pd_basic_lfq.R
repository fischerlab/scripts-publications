#####
# IP script adapted from TMT
# Fischer lab
# Dependencies: Fischer lab helper functions
# This script is for basic extraction of expression values and limma statistical analysis from PD proteins files

# Define variables
### load PD data
D.exp <- read.delim("~/", header=TRUE, sep="\t", as.is=TRUE)
dim(D.exp) 

#number of minimal spectral counts
sc <- 2
#dataset prefix
pre <- "file-name"
# use intensities (1) or spectral counts (2) as weights for limma
w <- 1
# type of input file (1) PD21, (2) MQ15, (3) Gygi core quant (4) other
qfiletype <- 4



# load helper functions
source("~/Dropbox (Partners HealthCare)/R/R-helper/map2col.R")
source("~/Dropbox (Partners HealthCare)/R/R-helper/plotPairwiseCorr.R")
source("~/Dropbox (Partners HealthCare)/R/R-helper/mypairs.R")
library("limma")
library("dplyr")



# Remove protein contaminants
sum(D.exp$Contaminant=="True")
D.exp <- D.exp[!D.exp$Contaminant=="True",]

# Identify  reporter ion channels
colnames(D.exp)
if(qfiletype==1){
  channels <- grep("^Abundance.F..(126|127|128|129|130|131)", colnames(D.exp), value=TRUE)
} else if(qfiletype==2){
  channels <- grep("^Reporter.intensity.[0-9].", colnames(D.exp), value=TRUE)
} else if(qfiletype==3){
  channels <- grep("*(126|127|128|129|130|131)", colnames(D.exp), value=TRUE)
} else{
  channels <- grep("^Abundance.F...", colnames(D.exp), value=TRUE)
}

# identify number of multiplexing channels
tmt <- length(channels)
tmt


# Identify contrasts
channels
temp <- vector(mode="character", length = length(channels))
for (i in 1:length(channels)) {
  temp[i] <- substr(channels[i],
                    regexpr("Sample.", channels[i])[1] + nchar("Sample."), # character index of cpd name start
                    nchar(channels[i])) # character index of cpd name end
}
temp2 <- unique(temp)
contrasts <- temp2[!grepl("DMSO", temp2)]
rm(i, temp, temp2)
contrasts



# Subset data frame
colnames(D.exp)
sel1 <- colnames(D.exp)[colnames(D.exp)%in% "Accession"]
sel2 <- grep("^Gene.Symbol", colnames(D.exp), value = TRUE)
sel3 <- grep("^Description", colnames(D.exp), value = TRUE)
sel4 <- grep("^Number.of.Unique.Peptides", colnames(D.exp), value = TRUE)
sel5 <- c(sel1, sel2, sel3, sel4, channels)
sel5
D <- D.exp[, sel5]
rownames(D) <- D[,1]
rm(sel1,sel2,sel3,sel4,sel5)

colnames(D)


###### cleanup data ######
# Remove data with unique peptides below threshold
sum(!D$Number.of.Unique.Peptides >= sc)
D <- D[D$Number.of.Unique.Peptides >= sc,]

colnames(D)


# add imputation method here if needed


######  Normalize and scale data. Create relative abundance matrix #############
# Normalize and scale data to create relative abundance matrix
x <- colSums(D[,c(5:(5+tmt-1))])
x <- x/max(x)
y <- t(apply(D[,c(5:(5+tmt-1))],x, "/", x))
Dn <- D
Dn[,c(5:(5+tmt-1))] <- y


# create matrix for abundance and plot unscaled data
A <- as.matrix(Dn[,c(5:(5+tmt-1))])
summary(A)
png(filename = paste("plots/",pre,"_unscaled_intensities_pairs.png", sep = ""), width = 1200, height = 1200)
mypairs(log10(A))
dev.off()


# Scaling protein TMT channels to 100
x <- rowSums(A)
ZZ <- A/x
zz <- ZZ*100
# create matrix for relative abundance
a <- zz
rm(zz); rm(ZZ); rm(x)


######## Use limma to generate moderated t-statistic p-values ######
### Define experimental factors
colnames(a)
contrast <- "Sample"
#a <- a_backup
colnames(a)
a <- a[,c(1:6)]
treat <- factor(ifelse(grepl(contrast, colnames(a)),"DRUG","CTRL"), levels=c("CTRL","DRUG"))
ctrl <- factor(ifelse(grepl("DMSO",colnames(a)),"DMSO","TREAT"), levels=c("DMSO", "TREAT"))
treat
ctrl

######## LIMMA MODERATED T-STATISTICS #####
#design <- model.matrix(~ treat + ctrl) # define matrix
design <- model.matrix(~ treat) # define matrix
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



######Add gene names to tt for PD1 data - adapt for Gygi/MQ#######
combineBy <- "Accession"
identifier <- D[,c(1,2)]
tt2 <- tt
colnames(tt2) <- c("Accession", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
ids <- unique(c(tt2[,combineBy], identifier[,combineBy]))
tt3 <- cbind(tt2[match(ids, tt2[,combineBy]),], identifier[match(ids, identifier[,combineBy]),])
tt[,1] <- tt3[,9]
rm(tt2);rm(tt3)

#Write hit table
write.csv(tt, file=paste("tables/",pre, "_Limma_output_",contrast,"_vs_CTRLs.csv"))


########################################################################
###### MAKE PLOTS FOR LIMMA ANALYSIS ###################################
########################################################################
# Plot & Hits by p.value
# Define p-value cut-off
pval <- 0.01
lab_idx <- log10(pval)+0.4
# Define minimal logFC
lfc <- 1.25 # 100% change
lfc_idx <- log2(lfc) - 0.2
# find y-axis value
maxFC <- max(abs(tt$logFC))
#maxFC <- 2.5

pdf(file = paste("plots/",pre,"_hits_",contrast,"_Pvalue.pdf"), width = 12, height = 12)
sel <- tt$logFC <= -log2(lfc) & tt$P.Value <= pval | tt$logFC >= log2(lfc) & tt$P.Value <= pval
summary(sel)
#sel2 <- tt$logFC <= -log2(lfc) & tt$P.Value <= pval & tt$P.Value >= pval
#summary(sel2)
#plot
par(mar=c(5,4,4-3,2)+.1, cex=1.8)
plot(log10(tt$P.Value),tt$logFC,  xlab="-log10 P value", ylab=paste("logFC (log2  ", contrast," / DMSO)"),
     ylim=c(-(maxFC),maxFC), 
     cex=0.8,
     pch=20, col=ifelse(sel, "red", "#22222222"))
text(log10(tt[sel,"P.Value"]), tt[sel,"logFC"],adj = c(-0.25,0.55), labels=tt[sel,"ID"], col="black", cex=1.0)
#text(log10(tt[sel2,"adj.P.Val"]), tt[sel2,"logFC"],adj = c(-0.25,0.55), labels=tt[sel2,"ID"], col="black", cex=0.8)
abline(v=log10(pval), lty=2)       
abline(h=c(-1,1)*log2(lfc), lty=2) 
legend(pch=20, "topleft", bty="n", legend=paste("Hits: p < ", pval, "; logFC > ",lfc, "-fold"), cex=1.0, col="red")
text(-3.75, 0.2, labels=paste(lfc, "-fold up-regulated  "), cex = 1)
text(-3.72,-0.2, labels=paste(lfc, "-fold down-regulated"), cex = 1)
text(-0.6,-0.8, labels = paste("p = ",pval), cex = 1)
dev.off()



