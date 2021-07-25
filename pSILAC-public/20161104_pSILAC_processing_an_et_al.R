####
# Process time-course mass spec experiment
# pSILAC with T6=6h, T8=10h and T16=16h time-points.
####

library("RColorBrewer")
library("gplots")
library("data.table")
library("Biobase")
library("limma")
source("~/R-helper/map2col.R")
source("~/R-helper/selbach_Thalf_func.R")


#### globals
inFile <- "rawdata/20141105_T6_T8_T16_LFQ_proteinGroups.txt"

### load raw data
D <- read.delim(inFile, header=TRUE, sep="\t", as.is=TRUE)

### remove false positives
isFalsePositive <- D$Only.identified.by.site=="+" | D$Reverse=="+" | D$Potential.contaminant=="+" #| D$Ratio.H.L.count<2
summary(isFalsePositive)

#   Mode   FALSE    TRUE    NA's
#logical    6328     528       0 # "rawdata/20141105_T6_T8_T16_LFQ_proteinGroups.txt" #Orbitrap experiment
D <- D[!isFalsePositive,]

##### define true positive - validated targets #####
isTruePositive <- D$Gene.names %in% c("CSNK1A1;CSNK1A1L", "ZFP91;ZFP91-CNTF", "RNF166", "ZNF692") # CSNK1A1, ZFP91, RNF166, ZNF692
sum(isTruePositive) # 2
#

#### flip label swap experiments and bind H/L, counts, int, LFQint data as matrix
# Labels:
# C = count of peptides
# M = H/L ratios
# A = Intensities
# Al??= LFQ ratios Light peptides
# Ah = LFQ ratios Heavy peptides
# Ai = sum of LFQ H+L ratios
# ...just T6 
#counts
nms <- grep("^Ratio.H.L.count.T6",colnames(D),value=TRUE)
C6 <- as.matrix(D[,nms])
#H/L ratios
nms <- grep("^Ratio.H.L.T6",colnames(D),value=TRUE)
M6 <- log2(as.matrix(D[,nms]))
if(any(i <- grepl("_LS$",nms)))
  M6[,i] <- -1 *M6[,i]
# intensities
nms <- grep("^Intensity.T6",colnames(D),value=TRUE)
A6 <- log10(D[,nms]+1)
# LFQ intensities
nms <- grep("^LFQ.intensity.L.T6",colnames(D),value=TRUE)
Al6tmp <- log10(as.matrix(D[,nms]+1))
Al6 <- Al6tmp
nms <- grep("^LFQ.intensity.H.T6",colnames(D),value=TRUE)
Ah6tmp <- log10(as.matrix(D[,nms]+1))
Ah6 <- Ah6tmp
# exchange LS experiments for LFQ quants
colnames(Al6tmp)
nms <- grep("^LFQ.intensity.L.T6",colnames(Al6tmp),value=TRUE)
if(any(i <- grepl("_LS$",nms)))
  Al6[,i] <- Ah6tmp[,i]
if(any(i <- grepl("_LS$",nms)))
  Ah6[,i] <- Al6tmp[,i]
Ai6 <- Al6 + Ah6

# ...just T8
nms <- grep("^Ratio.H.L.count.T8",colnames(D),value=TRUE)
C8 <- as.matrix(D[,nms])
nms <- grep("^Ratio.H.L.T8",colnames(D),value=TRUE)
M8 <- log2(as.matrix(D[,nms]))
if(any(i <- grepl("_LS$",nms)))
  M8[,i] <- -1 *M8[,i]
nms <- grep("^Intensity.T8",colnames(D),value=TRUE)
A8 <- log10(as.matrix(D[,nms]+1))
# LFQ intensities
nms <- grep("^LFQ.intensity.L.T8",colnames(D),value=TRUE)
Al8tmp <- log10(as.matrix(D[,nms]+1))
Al8 <- Al8tmp
nms <- grep("^LFQ.intensity.H.T8",colnames(D),value=TRUE)
Ah8tmp <- log10(as.matrix(D[,nms]+1))
Ah8 <- Ah8tmp
# exchange LS experiments for LFQ quants
colnames(Al8tmp)
nms <- grep("^LFQ.intensity.L.T8",colnames(Al8tmp),value=TRUE)
if(any(i <- grepl("_LS$",nms)))
  Al8[,i] <- Ah8tmp[,i]
if(any(i <- grepl("_LS$",nms)))
  Ah8[,i] <- Al8tmp[,i]
Ai8 <- Al8 + Ah8

# ...just T16
nms <- grep("^Ratio.H.L.count.T16",colnames(D),value=TRUE)
C16 <- as.matrix(D[,nms])
nms <- grep("^Ratio.H.L.T16",colnames(D),value=TRUE)
M16 <- log2(as.matrix(D[,nms]))
if(any(i <- grepl("_LS$",nms)))
  M16[,i] <- -1 *M16[,i]
nms <- grep("^Intensity.T16",colnames(D),value=TRUE)
A16 <- log10(as.matrix(D[,nms]+1))
# LFQ intensities
nms <- grep("^LFQ.intensity.L.T16",colnames(D),value=TRUE)
Al16tmp <- log10(as.matrix(D[,nms]+1))
Al16 <- Al16tmp
nms <- grep("^LFQ.intensity.H.T16",colnames(D),value=TRUE)
Ah16tmp <- log10(as.matrix(D[,nms]+1))
Ah16 <- Ah16tmp
# exchange LS experiments for LFQ quants
colnames(Al6tmp)
nms <- grep("^LFQ.intensity.L.T16",colnames(Al16tmp),value=TRUE)
if(any(i <- grepl("_LS$",nms)))
  Al16[,i] <- Ah16tmp[,i]
if(any(i <- grepl("_LS$",nms)))
  Ah16[,i] <- Al16tmp[,i]
Ai16 <- Al16 + Ah16

# ...all
#counts
C <- cbind(C6,C8,C16)
#ratios
M <- cbind(M6,M8,M16)
#intensities
A <- cbind(A6,A8,A16)
#LFQ intensities
Al <- cbind(Al6,Al8,Al16)
Ah <- cbind(Ah6,Ah8,Ah16)
Ai <- cbind(Ai6,Ai8,Ai16)

######### LOOK FOR TRANSLATIONAL CHANGES #######
sel1 <- rowSums(C[,9:12]>=2) == ncol(C[,9:12]) #| rownames(Mn) %in% gnid.truePositives #| xx$SDdmso <= 0.5 & xx$SDlena <= 0.5
summary(sel1)

# remove NA
colnames(C[,9:12])
rownames(Ah) <- rownames(Al) <- D$Gene.names

ZZ <- cbind(Ah[,9:12],Al[,9:12],C)
ZZ <- ZZ[sel1,]
summary(ZZ)
colnames(ZZ)
ZZ[ZZ==0] <- NA
ZZ <- na.omit(ZZ)

#create matrix for T16 Heavy proteins
Z <- ZZ[,1:4]
colnames(Z)
summary(Z)
hist(Z[,1])
tab <- cbind(rowMeans(Z[,c(1:2)]),
             rowMeans(Z[,c(3:4)]))
temp1 <- as.matrix(tab[,1] - tab[,2])
hist(temp1, breaks=100)
### LOOK AT LIGHHT T16 proteins ####
Y <- ZZ[,5:8]
colnames(Y)
summary(Y)
hist(Y[,1])
tab <- cbind(rowMeans(Y[,c(1:2)]),
             rowMeans(Y[,c(3:4)]))
temp2 <- as.matrix(tab[,1] - tab[,2])
hist(temp2, breaks=100)

plot(temp1, temp2)
temp3 <- as.matrix(temp1 - temp2)

### extract experimental factors
treat <- factor(ifelse(grepl("DMSO",colnames(M)),"DMSO","LENA"), levels=c("DMSO","LENA"))
hours <- factor(sub("^.*(T[0-9]+).*$","\\1",colnames(M)), levels=c("T6","T10","T16"))
hours <- factor(c("T6","T6","T6","T6","T6","T6","T10", "T10","T16","T16","T16","T16"), levels=c("T6","T10","T16"))
hoursNum <- as.numeric(sub("^.*T([0-9]+).*$","\\1",colnames(M)))
hoursNum <- c(6,6,6,6,6,6,10,10,16,16,16,16)

####### Compare H/L distribution before normalization#######
dL <- lapply(seq.int(ncol(M)), function(i) density(M[,i],na.rm=TRUE,bw=0.1))
xlims <- c(-4,4)
ylims <- c(0,max(sapply(dL, function(d) max(d$y))))
cols <- brewer.pal(3,'Set1')[hours]
ltys <- c(3,1)[treat]
pdf("plots/20160204_EDXFig1_RatioHL_distributions_beforeNorm.pdf", pointsize=16)
plot(0:1, 0:1, type="n", xlab="Ratio log2(H/L)", ylab="Density", main="", xlim=xlims, ylim=ylims) 
for(i in c(2:3,5:12)) # remove T6_LS (label swap)
  lines(dL[[i]], col=cols[i], lty=ltys[i], lwd=2)
legend(x="topright",bty="n",lty=c(3,1),col=brewer.pal(3,'Set1')[rep(1:nlevels(hours),each=2)],
       legend=sprintf("%s %s",rep(levels(hours),each=2),levels(treat)))
dev.off()

###### NORMALIZE DATA - CLEAR OFFSET FOR LABEL SWAP EXPERIMENTS #####
### Requires normalization for each timepoint seperately ###
M6n <- normalizeBetweenArrays(M6)
M8n <- normalizeBetweenArrays(M8)
M16n <- normalizeBetweenArrays(M16)
Mn <- cbind(M6n,M8n,M16n)
rownames(Mn) <- D$Gene.names

### check normalized data "Mn"
dL <- lapply(seq.int(ncol(Mn)), function(i) density(Mn[,i],na.rm=TRUE,bw=0.1))
xlims <- c(-4,4)
ylims <- c(0,max(sapply(dL, function(d) max(d$y))))
cols <- brewer.pal(3,'Set1')[hours]
ltys <- c(3,1)[treat]
pdf("plots/20160204_EDXFig1_RatioHL_distributions_afterNorm.pdf", pointsize=16)
plot(0:1, 0:1, type="n", xlab="Ratio log2(H/L)", ylab="Density", main="", xlim=xlims, ylim=ylims) 
for(i in c(2:3,5:12)) # remove T6_LS (label swap)
  lines(dL[[i]], col=cols[i], lty=ltys[i], lwd=2)
legend(x="topright",bty="n",lty=c(3,1),col=brewer.pal(3,'Set1')[rep(1:nlevels(hours),each=2)],
       legend=sprintf("%s %s",rep(levels(hours),each=2),levels(treat)))
dev.off()
hist(Mn[,10],breaks=200, xlim = c(-5,5))
summary(Mn)
##################################



############### Average replicates and look at delta Ratios ###########
sel <- apply(Mn,1,function(x) all(is.finite(x))); sum(sel) # 2837
#average replicates T16
x16 <- rowMeans(subset(Mn, select = c("Ratio.H.L.T16_DMSO1", "Ratio.H.L.T16_DMSO2")), na.rm = TRUE)
y16 <- rowMeans(subset(Mn, select = c("Ratio.H.L.T16_LENA2", "Ratio.H.L.T16.LENA1")), na.rm = TRUE)

#explorative plot
png("plots/20160305_T16_deltaHL_average.png", height=1200, width=1200, pointsize=25)
plot(x16[sel], y16[sel], xlab = "T16 H/L DMSO", ylab = "T16 H/L LENA", pch="*", cex=0.8)#, col="#22222222")
abline(a=0,b=1,lty=3)
abline(a=0.66,b=1,lty=3)
abline(a=-0.66,b=1,lty=3)
#highlight True positives
isTruePositive2 <- D$Gene.names %in% c("CSNK1A1;CSNK1A1L", "ZFP91;ZFP91-CNTF")
isTruePositive3 <- D$Gene.names %in% c("CSNK1A1;CSNK1A1L")
points(x16[isTruePositive2], y16[isTruePositive2], pch="*", cex=1.2, col=82)
points(x16[isTruePositive3], y16[isTruePositive3], pch="*", cex=1.2, col="blue")

legend(x="bottomright", bty="n", pch="*", cex=1.2, legend="CSNK1A1", col=82)
dev.off()

test <- as.matrix(cbind(x16[sel], y16[sel]))
sd(test)


##############################################
#
#Fit half life according to Selbach et al.
#
################

### fit decay
# expect to see greater ratios at later compared to earlier (more turnover at later timepoint)
sel <- apply(Mn,1,function(x) all(is.finite(x))); sum(sel) # 2837
table(isTruePositive,sel) # lost one of the two due to NaN values
colnames(Mn)

# visualize time dependence # FIgure 1b
png("plots/20160315_RatioHL_scatter_T6-vs-later.png", height = 1200, width = 1200, pointsize=20)
par(mar=c(5,4,4-3,2)+.1)
plot(rowMeans(M6n[sel,]), rowMeans(M8n[sel,]), pch="*", xlab="T6 Ratio log2(H/L)", ylab="TX Ratio log2(H/L)",
     xlim = c(-2.8,1.8),
     ylim = c(-1.8,3))
points(rowMeans(M6n[sel,]), rowMeans(M16n[sel,]), pch="*", col="darkred")
abline(a=0,b=1)
legend(x="bottomright", bty="n", pch="*", col=c("black","darkred"), legend=c("T8","T16"))
#points(rowMeans(M6n[sel & isTruePositive,,drop=FALSE]), rowMeans(M8n[sel & isTruePositive,,drop=FALSE]), pch=1)
#points(rowMeans(M6n[sel & isTruePositive,,drop=FALSE]), rowMeans(M16n[sel & isTruePositive,,drop=FALSE]), pch=1)
#legend(x="topleft", bty="n", pch=1, legend=D$Gene.names[sel & isTruePositive])
dev.off()



####### LOOK AT HISTOGRAMS OF RATIIOS Figure 1c
colnames(Mn)
Mh6 <- rowMeans(Mn[,c(2,3,5,6)])
Mh10 <- rowMeans(Mn[,c(7,8)])
Mh16 <- rowMeans(Mn[,c(9:12)])

png("plots/20160419_RatioHL_Hist_T6-vs-later.png", height = 1200, width = 1200, pointsize=30)
hist(Mh6, col=rgb(0,0,1,1/4), xlim=c(-3,3), ylim = c(0,500), breaks=50, cex.lab=1.5, cex.axis=1.5)
hist(Mh10, col=rgb(0,1,0,1/4), xlim=c(-3,3), ylim = c(0,500), breaks=50, add=T)
hist(Mh16, col=rgb(1,0,0,1/4), xlim=c(-3,3), ylim = c(0,500), breaks=50, add=T)
legend(x="topright", bty="n", pch=15, col=c(rgb(0,0,1,1),rgb(0,1,0,1),rgb(1,0,0,1)), legend=c("T6","T10","T16"), cex=1.5)
dev.off()

############


Thalf <- apply(Mn,1,calcThalf,ti=hoursNum,tcc=100)
hist(Thalf,100)

#### print log plot of some examplary genes #####
plot(hoursNum, Mn[D$Gene.names %in% c("UBA1")],
     ylim = c(-2.2,4))
points(hoursNum, Mn[D$Gene.names %in% c("DDB1")], col="blue")
#points(hoursNum, Mn[D$Gene.names %in% c("CSNK1A1;CSNK1A1L")], col="blue")
points(hoursNum, Mn[D$Gene.names %in% c("ZFP91;ZFP91-CNTF")], col="red")
points(hoursNum, Mn[D$Gene.names %in% c("RRM2")], pch=3, col="red")
points(hoursNum, Mn[D$Gene.names %in% c("HIST1H1C")], pch=3, col="red")
points(hoursNum, Mn[D$Gene.names %in% c("VPRBP")], pch=3, col="red")

Thalf.DMSO <- apply(Mn[,treat=="DMSO"],1,calcThalf,ti=hoursNum[treat=="DMSO"],tcc=100)
Thalf.LENA <- apply(Mn[,treat=="LENA"],1,calcThalf,ti=hoursNum[treat=="LENA"],tcc=100)
hist(Thalf.DMSO, 100, xlab="Protein half-life (hours)", main="", col="gray")
hist(Thalf.LENA, 100, xlab="Protein half-life (hours)", main="", col="gray")


pdf("plots/20160110_Thalf_T6-T8-T16_histogram.pdf", height=6, width=8, pointsize=16)
par(mar=c(5,4,4-3,2)+.1)
hist(Thalf, 100, xlab="Protein half-life (hours)", main="", col="gray")
abline(v=median(Thalf,na.rm=TRUE), lty=2)
text(x=median(Thalf,na.rm=TRUE)-par('cxy')[1], y=par('usr')[4], adj=c(1,1), labels=sprintf("%.1f h",median(Thalf,na.rm=TRUE)), xpd=NA)
rug(Thalf[isTruePositive],col="red")
legend(x="topright", text.col="red", legend="true positive", bty="n")
dev.off()

png("plots/20141106_Thalf_T6-T8-T16_DMSO_vs_LENA_scatter.png", height=800, width=1600, pointsize=20)
par(mar=c(5,4,4-3,2)+.1)
z <- rowMeans(A,na.rm=T)
sel1 <- rowSums(C>=1) == ncol(C); summary(sel1) # at least 1 peptide  per sample
sel2 <- rowSums(C>=2) == ncol(C); summary(sel2) # at least 2 peptides per sample
#table(isTruePositive, sel1)
#table(isTruePositive, sel2)
isTruePositive2 <- D$Gene.names %in% c("CSNK1A1;CSNK1A1L", "ZFP91", "ZFP91;ZFP91-CNTF")
par(mar=c(5,4,4-3,2)+.1, mfrow=c(1,2))
plot(Thalf.DMSO[sel1], Thalf.LENA[sel1], pch=21, cex=0.8, col="black", bg=map2col(z)[sel1],
     xlab="Protein half-life (DMSO, hours)", ylab="Protein half-life (LENA, hours)")
points(Thalf.DMSO[isTruePositive2], Thalf.LENA[isTruePositive2], pch="*", cex=1.2)
legend(x="topleft", bty="n", legend=sprintf(">=1 peptides per sample\n(n=%d)",sum(sel1)))
legend(x="bottomright", bty="n", pch="*", cex=1.2, legend="CSNK1A1")
abline(a=0,b=1,lty=3)
plot(Thalf.DMSO[sel2], Thalf.LENA[sel2], pch=21, cex=0.8, col="black", bg=map2col(z)[sel2],
     xlab="Protein half-life (DMSO, hours)", ylab="Protein half-life (LENA, hours)")
points(Thalf.DMSO[isTruePositive2], Thalf.LENA[isTruePositive2], pch="*", cex=1.2)
legend(x="topleft", bty="n", legend=sprintf(">=2 peptides per sample\n(n=%d)",sum(sel2)))
legend(x="bottomright", bty="n", pch="*", cex=1.2, legend="CSNK1A1")
abline(a=0,b=1,lty=3)
dev.off()
#Quality cutoff for 1 peptide/sample is sufficient.

### calculated delta H/L ratios (len - dmso)
MM  <- Mn[,treat=="LENA"] - Mn[,treat=="DMSO"]
sel2 <- rowSums(C>=2) == ncol(C); summary(sel2) # at least 2 peptides per sample


png("plots/20141106_changeOfRatio_T16_rep1-vs-repl2_scatter.png", height=800, width=800, pointsize=20)
par(mar=c(5,4,4-3,2)+.1)
z <- rowMeans(A,na.rm=T)
sel <- rowSums(C>=2) == ncol(C); sum(sel) # at least 2 peptides per sample # 2230
#sel <- rowSums(C>=1) == ncol(C); sum(sel) # at least 1 peptide per sample # 2837
x <- rowMeans(MM[sel,"Ratio.H.L.T16.LENA1",drop=FALSE])
y <- rowMeans(MM[sel,"Ratio.H.L.T16_LENA2",drop=FALSE])
plot(x, y, pch=21, cex=0.8, col="black", bg=map2col(z)[sel],
     xlab="Change of H/L ratio (T16 replicate 1)", ylab="Change of H/L ratio (T16 replicate 2)")
xtp <- rowMeans(MM[isTruePositive,"Ratio.H.L.T16.LENA1",drop=FALSE],na.rm=TRUE)
ytp <- rowMeans(MM[isTruePositive,"Ratio.H.L.T16_LENA2",drop=FALSE],na.rm=TRUE)
points(xtp, ytp, pch="*", cex=1.2)
legend(x="topleft", bty="n", legend=sprintf("R = %.3f",cor(x,y,use="pair")))
legend(x="bottomright", bty="n", pch="*", cex=1.2, legend="true positive")
dev.off()


# ... compare to half-time plot
sel8 <- x < -0.2 & y < -0.2; sum(sel8)
sel9 <- x > 0.2 & y > 0.2; sum(sel9)

sel2 <- rowSums(C>=2) == ncol(C); summary(sel2) # at least 2 peptides per sample
plot(Thalf.DMSO[sel2], Thalf.LENA[sel2], pch=21, cex=0.8, col="black", bg=map2col(z)[sel2],
     xlab="Protein half-life (DMSO, hours)", ylab="Protein half-life (LENA, hours)")
legend(x="topleft", bty="n", legend=sprintf(">=2 peptides per sample\n(n=%d)",sum(sel2)))
abline(a=0,b=1,lty=3)
points(Thalf.DMSO[names(x)[sel8]], Thalf.LENA[names(x)[sel8]], pch="*", col="blue")
points(Thalf.DMSO[names(x)[sel9]], Thalf.LENA[names(x)[sel9]], pch="*", col="green")

#D[names(x)[sel8],c("Protein.IDs","Gene.names","Fasta.headers")]
D[names(x)[sel8],c("Protein.IDs","Gene.names")]
D[names(x)[sel9],c("Protein.IDs","Gene.names")]

selGene <- names(x)[sel8][1]
selGene <- names(x)[sel8][5]

plot(hoursNum, M[selGene,], col=brewer.pal(3,'Set1')[treat], xlab="Measurement time (hours)", ylab="Ratio log2(H/L)", log="")

pdf("plots/20160226_time-vs-ratio_linearFits_topInhibitedProteinsAnTruePositives.pdf", height=9, width=12, pointsize=16)
par(mfrow=c(3,4), mar=c(5,4,4-3,2)+.1)
for(i in 1:(12-sum(isTruePositive))) {
  selGene <- names(x)[sel8][i]
  plot(hoursNum, M[selGene,], col=brewer.pal(3,'Set1')[treat], xlab="Measurement time (hours)", ylab="Ratio log2(H/L)", log="")
  abline(f1 <- lm(M[selGene,] ~ hoursNum, subset=treat=="DMSO"), col=brewer.pal(3,'Set1')[1])
  abline(f2 <- lm(M[selGene,] ~ hoursNum, subset=treat=="LENA"), col=brewer.pal(3,'Set1')[2])
  legend(x="topleft", bty="n", legend=D[selGene,'Gene.names'])
  legend(x="bottomright", bty="n", text.col=brewer.pal(3,'Set1')[1:2], legend=sprintf("%s (%.2f)",levels(treat),c(summary(f1)$adj.r.squared, summary(f2)$adj.r.squared)))
}
# true positives
for(selGene in rownames(M)[isTruePositive]) {
  plot(hoursNum, M[selGene,], col=brewer.pal(3,'Set1')[treat], xlab="Measurement time (hours)", ylab="Ratio log2(H/L)", log="")
  abline(f1 <- lm(M[selGene,] ~ hoursNum, subset=treat=="DMSO"), col=brewer.pal(3,'Set1')[1])
  abline(f2 <- lm(M[selGene,] ~ hoursNum, subset=treat=="LENA"), col=brewer.pal(3,'Set1')[2])
  legend(x="topleft", bty="n", legend=D[selGene,'Gene.names'])
  legend(x="bottomright", bty="n", text.col=brewer.pal(3,'Set1')[1:2], legend=sprintf("%s (%.2f)",levels(treat),c(summary(f1)$adj.r.squared, summary(f2)$adj.r.squared)))
}
dev.off()

########
isControls <- D$Gene.names %in% c("UBA1", "DDB1", "COPS5", "GAPDH") # CSNK1A1, ZFP91, RNF166, ZNF692
sum(isControls) # 4
colnames(M)
Mnls <- M[,c(2,3,5,6,7,8,9,10,11,12)]
Mnls[isControls]
hoursNum <- c(6,6,6,6,10,10,16,16,16,16)
treat <- factor(c("DMSO", "DMSO", "LENA", "LENA", "DMSO", "LENA", "DMSO", "DMSO", "LENA", "LENA"))
levels(treat) <- c("DMSO", "LENA")
treat


####### Sort data for r2 value of linear regression H/L vs. time ##### THIS IS THE FINAL QUALITY CRITERIUM USED - processing for figures is below
colnames(Mn)
#rownames(Mn) <- D$Gene.names
summary(Mn)
TempLm1 <- Mn[,c(2,3,5,6,7,8,9,10,11,12)]
colnames(TempLm1)
rownames(TempLm1)
sel <- rownames(TempLm1) == "ZFP91;ZFP91-CNTF"
summary(TempLm1)
TempLm1[sel]
sel <- apply(TempLm1,1,function(x) all(is.finite(x))) | rownames(TempLm1) %in% "ZFP91;ZFP91-CNTF"
###
#colnames(TempLm1)
#n <- apply(as.matrix(TempLm1), 1, function(x) length(is.na(x)[is.na(x)==FALSE]))
#TempLm1 <- as.matrix(TempLm1[n>7,])
####
TempLm1 <- as.matrix(TempLm1[sel,])
dim(TempLm1)
TempLm <- as.matrix(cbind(TempLm1, Rdmso=rnorm(nrow(TempLm1)), Rlena=rnorm(nrow(TempLm1))))# , ID=rownames(TempLm1)))
#### impute missing values for ZFP91
selTemp <- rownames(TempLm) %in% "ZFP91;ZFP91-CNTF"
summary(selTemp)
selTemp
TempLm[2272,]
TempLm[2272,3] <- -1.0300562 
TempLm[2272,8] <- 0.5028062
TempLm[2272,10] <- 1.7078167
TempLm[selTemp]
colnames(TempLm)
rownames(TempLm)
hoursNum
treat
hoursNum <- c(6,6,6,6,10,10,16,16,16,16)
treat <- factor(c("DMSO", "DMSO", "LENA", "LENA", "DMSO", "LENA", "DMSO", "DMSO", "LENA", "LENA"), levels = c("DMSO", "LENA"))
treat
########### Calculate half-lives only if R2 of linear regression >= 0.9 as quality criterium. And ratio at all time points
# calculate r2 values and add to table
for(i in 1:nrow(TempLm1)){
  ftmpD <- lm(TempLm1[i,] ~ hoursNum, subset=treat=="DMSO")
  ftmpL <- lm(TempLm1[i,] ~ hoursNum, subset=treat=="LENA")
  TempLm[i,11] <- as.numeric(summary(ftmpD)$adj.r.squared)
  TempLm[i,12] <- summary(ftmpL)$adj.r.squared
}
## select data points that fullfil r2 criteria
selLm <- TempLm[,11] >= 0.9 & TempLm[,12] >= 0.9
summary(selLm)
# calculate Thalf
colnames(TempLm)
Thalf.DMSO <- apply(TempLm[,c(1,2,5,7,8)],1,calcThalf,ti=c(6,6,10,16,16),tcc=150)
Thalf.LENA <- apply(TempLm[,c(3,4,6,9,10)],1,calcThalf,ti=c(6,6,10,16,16),tcc=150)


filtThalf <- data.frame(Thalf.DMSO, Thalf.LENA, rownames(TempLm))
selGenFit <- filtThalf$rownames.TempLm. %in% c("CSNK1A1;CSNK1A1L", "ZFP91;ZFP91-CNTF")
filtThalf$Thalf.DMSO[selGenFit]
filtThalf$Thalf.LENA[selGenFit]

####### Figure 1e ###########
plot(filtThalf$Thalf.DMSO[selLm], filtThalf$Thalf.LENA[selLm], pch=16, col="#22222222")
abline(a=0,b=1,lty=3)
#3 STDEV abline
abline(-(3*0.35),1, lty=3)
abline((3*0.35),1, lty=3)
#5 STDEV abline
abline(-(5*0.35),1, lty=3)
abline((5*0.35),1, lty=3)

summary(selGenFit)
points(filtThalf$Thalf.DMSO[selGenFit], filtThalf$Thalf.LENA[selGenFit], col="red", pch=16, cex=1.0)
text(filtThalf[selGenFit, "Thalf.DMSO"] + par('cxy')[1], filtThalf[selGenFit, "Thalf.LENA"], adj=c(0,0.5), labels=filtThalf[selGenFit, "rownames.TempLm."], col="black")
# Look at hits in lower cut-off
hitsUp <- (1*filtThalf$Thalf.DMSO+(5*0.35)) <= filtThalf$Thalf.LENA 
hitsDown <- (1*filtThalf$Thalf.DMSO-(5*0.35)) >= filtThalf$Thalf.LENA 

hitsUp <- (1*filtThalf$Thalf.DMSO+(3*0.35)) <= filtThalf$Thalf.LENA 
hitsDown <- (1*filtThalf$Thalf.DMSO-(3*0.35)) >= filtThalf$Thalf.LENA 
filtThalf$rownames.TempLm.[hitsUp & selLm]
filtThalf$rownames.TempLm.[hitsDown & selLm]
#points(filtThalf$Thalf.DMSO[hitsUp], filtThalf$Thalf.LENA[hitsUp], col="green", pch=1)
#points(filtThalf$Thalf.DMSO[hitsDown], filtThalf$Thalf.LENA[hitsDown], col="green", pch=1)
#text(filtThalf[hitsDown, "Thalf.DMSO"] + par('cxy')[1], filtThalf[hitsDown, "Thalf.LENA"], adj=c(0,0.5), labels=filtThalf[hitsDown, "rownames.TempLm."], col="black")
#text(filtThalf[hitsUp, "Thalf.DMSO"] + par('cxy')[1], filtThalf[hitsUp, "Thalf.LENA"], adj=c(0,0.5), labels=filtThalf[hitsUp, "rownames.TempLm."], col="black")
dev.off()
temp <-  filtThalf$Thalf.LENA[selLm] - filtThalf$Thalf.DMSO[selLm]
temp <- data.frame(temp, filtThalf$rownames.TempLm.[selLm])
####### Calculate SD and check for normal distribution
quantile(temp[,1], c(0.1587, 0.5, 0.8413), na.rm=TRUE)
sd(temp[,1])
mad(temp[,1])
#left SD = -0.212, right SD = 0.239 according to quartiles.
# conventional SD = 0.331
hist(temp[,1], breaks=50)
shapiro.test(temp[,1]) #W = 0.89393, p-value < 2.2e-16

##### END OF ## Calculate half-lives only if R2 of linear regression >= 0.9 as quality criterium. And ratio at all time points



####### Revised Figure 1e ###########
pdf("plots/20161114_EDFig_XX_newHalfLife_cut-offs.pdf", height=12, width=12, pointsize=16, useDingbats = FALSE)
plot(filtThalf$Thalf.DMSO[selLm], filtThalf$Thalf.LENA[selLm], pch=16, col="#22222222")
abline(a=0,b=1,lty=3)
#3 STDEV abline
abline(-(3*0.35),1, lty=3, col="blue")
abline((3*0.35),1, lty=3, col="blue")
#5 STDEV abline
abline(-(5*0.35),1, lty=3, col="red")
abline((5*0.35),1, lty=3, col="red")
summary(selGenFit)
points(filtThalf$Thalf.DMSO[selGenFit], filtThalf$Thalf.LENA[selGenFit], col="red", pch=16, cex=1.0)
text(filtThalf[selGenFit, "Thalf.DMSO"] + par('cxy')[1], filtThalf[selGenFit, "Thalf.LENA"], adj=c(0,0.5), labels=filtThalf[selGenFit, "rownames.TempLm."], col="black")
# Look at hits in lower cut-off
# 5 Stdev
hitsUp <- (1*filtThalf$Thalf.DMSO+(5*0.35)) <= filtThalf$Thalf.LENA 
hitsDown <- (1*filtThalf$Thalf.DMSO-(5*0.35)) >= filtThalf$Thalf.LENA 
# 3 Stdev
#hitsUp <- (1*filtThalf$Thalf.DMSO+(3*0.35)) <= filtThalf$Thalf.LENA 
#hitsDown <- (1*filtThalf$Thalf.DMSO-(3*0.35)) >= filtThalf$Thalf.LENA 

filtThalf$rownames.TempLm.[hitsUp & selLm]
filtThalf$rownames.TempLm.[hitsDown & selLm]

points(filtThalf$Thalf.DMSO[hitsUp & selLm], filtThalf$Thalf.LENA[hitsUp & selLm], col="red", pch=1)
points(filtThalf$Thalf.DMSO[hitsDown & selLm], filtThalf$Thalf.LENA[hitsDown & selLm], col="red", pch=1)
text(filtThalf[hitsDown & selLm, "Thalf.DMSO"] + par('cxy')[1], filtThalf[hitsDown & selLm, "Thalf.LENA"], adj=c(0,0.5), labels=filtThalf[hitsDown & selLm, "rownames.TempLm."], col="black")
text(filtThalf[hitsUp & selLm, "Thalf.DMSO"] + par('cxy')[1], filtThalf[hitsUp & selLm, "Thalf.LENA"], adj=c(1.5,0.5), labels=filtThalf[hitsUp & selLm, "rownames.TempLm."], col="black")
dev.off()
temp <-  filtThalf$Thalf.LENA[selLm] - filtThalf$Thalf.DMSO[selLm]
temp <- data.frame(temp, filtThalf$rownames.TempLm.[selLm])



######## FIGURE 1f #############
hist(temp[,1], breaks=50)
# calc deltaT16 with similar cut-offs
colnames(TempLm)
temp2 <- rowMeans(TempLm[,c(9,10)]) - rowMeans(TempLm[,c(7,8)])
temp2 <- temp2[selLm]

plot(temp[,1], temp2, pch=16, col="#22222222", xlim = c(-4,4), ylim = c(-1,1))
#abline(lm(temp2 ~ temp[,1]))
tempSel <- temp$filtThalf.rownames.TempLm..selLm. %in% c("CSNK1A1;CSNK1A1L", "ZFP91;ZFP91-CNTF")
points(temp[tempSel, "temp"], temp2[tempSel], pch=16, col="red")
text(temp[tempSel, "temp"] + par('cxy')[1], temp2[tempSel], adj=c(0,0.5), labels="CSNK1A1", col="black")
legend(x="topright", bty="n", legend=sprintf("R = %.3f",cor(temp[,1],temp2,use="pair",method="pearson")))
#sel <- temp$temp <= -1.5 & temp2 >= 0.5
#temp$filtThalf.rownames.TempLm..selLm.[sel]
#summary(sel)

###### END OF FIGURE 1f #############


########## FIGURE 1d ##########

pdf("plots/20160419_time-vs-ratio_linearFits_HitsAndCtrl.pdf", height=16.5, width=24, pointsize=25)
par(mfrow=c(2,3))
for(selGene in rownames(Mnls)[isTruePositive]) {
  plot(hoursNum, Mnls[selGene,], col=brewer.pal(3,'Set1')[treat], xlab="Measurement time (hours)", ylab="protein ratio log2(H/L)", log="",
       cex.axis=1.5, cex.lab=1.5, pch=1)
  abline(f1 <- lm(Mnls[selGene,] ~ hoursNum, subset=treat=="DMSO"), col=brewer.pal(3,'Set1')[1])
  abline(f2 <- lm(Mnls[selGene,] ~ hoursNum, subset=treat=="LENA"), col=brewer.pal(3,'Set1')[2])
  legend(x="topleft", bty="n", legend=D[selGene,'Gene.names'], cex=1.2)
  legend(x="bottomright", bty="n", text.col=brewer.pal(3,'Set1')[1:2], legend=sprintf("%s (%.2f)",levels(treat),c(summary(f1)$adj.r.squared, summary(f2)$adj.r.squared)), cex=1.2)
}
#controls
for(selGene in rownames(Mnls)[isControls]) {
  plot(hoursNum, Mnls[selGene,], col=brewer.pal(3,'Set1')[treat], xlab="Measurement time (hours)", ylab="Ratio log2(H/L)", log="",
       cex.axis=1.5, cex.lab=1.5, pch=1)
  abline(f1 <- lm(Mnls[selGene,] ~ hoursNum, subset=treat=="DMSO"), col=brewer.pal(3,'Set1')[1])
  abline(f2 <- lm(Mnls[selGene,] ~ hoursNum, subset=treat=="LENA"), col=brewer.pal(3,'Set1')[2])
  legend(x="topleft", bty="n", legend=D[selGene,'Gene.names'], cex=1.2)
  legend(x="bottomright", bty="n", text.col=brewer.pal(3,'Set1')[1:2], legend=sprintf("%s (%.2f)",levels(treat),c(summary(f1)$adj.r.squared, summary(f2)$adj.r.squared)), cex=1.2)
}
dev.off()
summary(f1)$adj.r.squared
########## END OF FIGURE 1d ##########





### visualize with fewer cut-offs (require finite numbers in each condition appears stringent enough)
#sel <- rowSums(head(is.finite(M))) >= 4
sel <- rowSums(C>=1) == ncol(C); sum(sel) # at least 1 peptide per sample
x <- Thalf.LENA - Thalf.DMSO; xlabel="Change of half-live (LENA-DMSO, hours)"
#y <- rowMeans(Ai[,c("LFQ.intensity.L.T16.LENA1","LFQ.intensity.L.T16.LENA2")] - Ai[,c("LFQ.intensity.L.T16_DMSO1","LFQ.intensity.L.T16_DMSO2")]); ylabel="Change in abundance (LENA-DMSO, log10)"
#y <- rowMeans(A[,c("Intensity.T16.LENA1","Intensity.T16_LENA2")] - Ai[,c("Intensity.T16_DMSO1","Intensity.T16_DMSO2")]); ylabel="Change in abundance (LENA-DMSO, log10)"
z <- rowMeans(MM[,c("Ratio.H.L.T16_LENA2","Ratio.H.L.T16.LENA1")],na.rm=TRUE); zlabel="Change of T16 Ratios (LENA-DMSO)"
sel2 <- sel & ((x < -1.25 & z >  0.4) | (x > 1.25 & z < -0.4)); sum(sel2)
sel3 <- sel & ((x < -1    & y < -0.2) | (x > 1.8  & y >  0.2)); sum(sel3)

png("plots/20160315_deltaHalftime-vs-delatRatio.png", height=1200, width=1200, pointsize=25)
#par(mar=c(5,4,4-3,2)+.1, mfrow=c(1,2))
plot(x[sel],z[sel],xlab=xlabel, ylab=zlabel, pch=20, col=map2col(rowMeans(A))[sel])
points(x[isTruePositive2], z[isTruePositive2], col="darkgreen", pch=1, cex=1.5)
#points(x[sel2],z[sel2],pch=1,cex=1)
#points(x[sel3],z[sel3],pch=3,cex=1)
legend(x="topright", bty="n", legend=sprintf("rho = %.3f",cor(x[sel],z[sel],use="pair",method="spearman")))
#legend(x="bottomleft", bty = "n", legend = c("H/L vs. half-life hits"
#                                             #, "Intensity vs. half-life hits"
#                                             ), cex = 0.8, pch = c(1, 3))
legend(x="bottomleft", bty = "n", legend = "CSNK1A1", cex = 0.8, pch =1, col = "darkgreen")
dev.off()

#plot(x[sel],y[sel],xlab=xlabel, ylab=ylabel, pch=20, col=map2col(rowMeans(A))[sel])
#points(x[sel2],y[sel2],pch=1,cex=1)
#points(x[sel3],y[sel3],pch="+",cex=1)
#points(x[isTruePositive], y[isTruePositive], pch=1, col="darkgreen", cex=1.5)
#legend(x="topright", bty="n", legend=sprintf("rho = %.3f",cor(x[sel],y[sel],use="pair",method="spearman")))
#legend(x="bottomright", bty = "n", legend = c("H/L vs. half-life hits", "Intensity vs. half-life hits"), cex = 0.8, pch = c(1, 3))
#dev.off()

listMulti <- data.frame(D[,c("Protein.IDs","Gene.names")], deltaHalflife=x, deltaAbundance=y, deltaT16ratio=z, sel.peptides=sel, sel.left=sel3, sel.right=sel2)
write.csv(listMulti, "output/20160111_deltaHalftime-vs-deltaAbundance-or-delatRatio_normalized.csv")

colnames(listMulti)
colnames(listMulti) <- c("Protein.IDs", "ID", "deltaHalflife", "deltaAbundance", "deltaT16ratio", "sel.peptides", "sel.left", "sel.right")

listAll <- merge(listLimma, listMulti, by = "ID")
write.csv(listAll, "output/20160111_hitlist_allTogether_Velos.csv")


