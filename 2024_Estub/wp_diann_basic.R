#####
# DIA-NN_wp script v0.0.1
# Fischer lab
# Dependencies: Fischer lab helper functions
# This script is for basic analysis of DIA-PASEF datasets from timsTOF processed with DIA-NN


#### Install DIA-NN R libraries ####
install.packages("devtools")
library(devtools)
install_github("https://github.com/vdemichev/diann-rpackage")

#Load packages
library("diann")
source("~/Dropbox (Partners HealthCare)/R/R-helper/mypairs.R")
library("gplots")
library("stringr")
library("limma")
library("dplyr")
library("readxl")


#Load data
df <- diann_load("")

w <- 1


#Preproces data
colnames(df)[1] <- "File.Name"
# Precursors x samples matrix filtered at 1% precursor and protein group FDR
precursors <- diann_matrix(df, pg.q = 0.01)
# Peptides without modifications - taking the maximum of the respective precursor quantities
peptides <- diann_matrix(df, id.header="Stripped.Sequence", pg.q = 0.01)
# Peptides without modifications - using the MaxLFQ algorithm
peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Stripped.Sequence", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
# Genes identified and quantified using proteotypic peptides
unique.genes <- diann_matrix(df, id.header="Genes", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01)
# Protein group quantities using MaxLFQ algorithm
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
dim(peptides.maxlfq)
dim(unique.genes)
dim(protein.groups)


### Explore data with maxLFQ quantified protein groups.
colnames(protein.groups)
D.exp <- as.data.frame(protein.groups)


#Define cut-off for min precuror 
sc <- 2

#Define samples and contrasts
samples <- c( "sample1", "sample2", "sample3", 
              "control1", "control2", "control3")  #
input.contrasts <- gsub("_rep[1|2|3]", "", samples)
input.channels <- input.contrasts  
colnames(D.exp) <- samples
head(D.exp)

# define prefix for output files
pre <- "add-experimental-details"
print(paste(sep="","filename prefix is ", pre))

###########################
#save workspace
save.image(file = "add-experimental-details.RData")  



# Identify number of runs
tmt <- length(colnames(D.exp))
tmt

print(paste(tmt, " runs identified", sep=""))

colnames(df)


D <- D.exp
colnames(D)

#Add back Meta Info
tmp <- as.data.frame(df[,c("Protein.Group", "Genes", "Q.Value", "PG.Q.Value")])
MetaInfo <- tmp[!duplicated(tmp$Protein.Group),]
head(MetaInfo)
occurences <- as.data.frame(table(unlist(tmp$Protein.Group)))
colnames(occurences) <- c("Protein.Group", "PrecursorForQuantTotal")
MetaInfo <- left_join(x = MetaInfo, y = occurences, by = "Protein.Group")
rm(tmp)
D <- cbind(rownames(D),D)
colnames(D) <- c("Protein.Group", colnames(D.exp))
tmp <- left_join(x = D, y = MetaInfo, by = "Protein.Group")
rownames(tmp) <- tmp$Protein.Group
colnames(tmp)
tmp <- cbind(tmp[,c(1,13,14,15,16,2:12)])
D <- tmp
rm(tmp)
colnames(D)

#Remove contaminants
contams<-read.delim(file = "~/Dropbox (Partners HealthCare)/R/R-helper/contams_list", header=T)
contams<-contams$Uniprot
D['Contaminant'] <- NA
D$Contaminant <- D$Protein.Group %in% contams | str_extract(D$Protein.Group,"[^;-]+") %in% contams
D <- D[!D$Contaminant,]
D <- subset(D, select = -c(Contaminant))

##explort all data before filtering
D.exp2 <- tibble::rownames_to_column(D.exp, "Accession.trimmed")
D.exp2 <- left_join(x = D.exp2, y = idmapping, by = "Accession.trimmed")

# Write hit table
write.csv(D.exp2, file=paste("./tables/",pre, "D-exp", ".csv"), row.names = F)


#Filter for minimum number of precursor
colnames(D)
sel <- D$PrecursorForQuantTotal >= sc*tmt
summary(sel)
D <- D[sel,]

#remove data with NAs or low sum of reporter ion intensities - threshold is experiment dependent..
sel <- as.matrix(rowSums(D[6:(6+tmt-1)])) >= #apply appropriate filters
summary(sel)
D <- subset(D, sel)
rm(sel)


# Normalize and scale data to create relative abundance matrix
x <- colSums(D[,c(6:(6+tmt-1))])
x <- x/max(x)
y <- t(apply(D[,c(6:(6+tmt-1))],x, "/", x))
Dn <- D
Dn[,c(6:(6+tmt-1))] <- y


#Skip in-house normalization - use if normalized in diann
x <- colSums(D[,c(6:(6+tmt-1))])
x <- x/max(x)
Dn <- D

# Plot normalization factors
png(filename = paste("./qc plots/",pre,"_TMT_normalization_factors.png", sep = ""), width = 800, height = 800)
bb <- barplot(x, ylab="relative abundance",
              ylim = c(0, 1.2*max(x)),
              names.arg = gsub("raw.PG.TMT.", "", colnames(D[,c(6:(6+tmt-1))])),
              cex.names = 1.2,
              cex.axis = 1.2,
              xlab = "TMT channels"
)
legend(pch=20, "topleft", bty="n", legend="Relative summed intensity for TMT channels / correction factors", cex=1.5, col="black")
text(x=bb,y=x, label=paste(round(x*100, digits=1),"%"), pos=3, cex=1.5, col="red")
dev.off()
rm(x); rm(y)

# Create matrix for abundance
A <- as.matrix(Dn[,c(6:(6+tmt-1))])

# Plot unscaled data
png(filename = paste("./qc plots/",pre,"_unscaled_intensities_pairs.png", sep = ""), width = 1200, height = 1200)
mypairs(log10(A))
dev.off()

# Scale matrix for abundance to 100
a <- log2((A/rowSums(A))*100) 

# Plot scaled data
png(filename = paste("./qc plots/",pre,"_scaled_intensities_pairs.png", sep = ""), width = 1200, height = 1200)
mypairs(a)
dev.off()


# Export scaled abundances
write.csv(file = paste("./scaled abundances/",pre,"_scaled_abundance_matrix.csv", sep=""), x=a, row.names = T)

#### Create a heatmap plot summarizing the data
#create color palette from green to red
my_palette <- colorRampPalette(c("royalblue1", "black", "sienna1"))(n=299)
### Heatmap for % RA
col_breaks = c(seq(0,10, length=100),               # for green
               seq(10.01, 20, length=100),            # for black
               seq(20.01, 30, length=100))             # for red
png(filename = paste("./qc plots/",pre,"_heatmap_all_channels.png", sep=""), width = 1800, height = 3600, res=300, pointsize = 8)
heatmap.2(2^a,
          na.rm = T,
          #cellnote = format(round(a, 2), nsmall = 2),
          main = "Relative abundance",
          notecol = "black",
          #density.info = "none",
          trace = "none",
          margins = c(12,9),
          col=my_palette,
          breaks = col_breaks,
          #hclustfun = function(x) hclust(a, method ="complete"),
          dendrogram = c("both"))
dev.off()
rm(my_palette,col_breaks)

#### Save image with sample clustering
Dmat = dist(t(a))
com.hclust = hclust(Dmat,method="complete")
png(filename = paste("./qc plots/",pre,"_sample_clustering.png", sep = ""), width = 1200, height = 1200)
plot(com.hclust,cex=1.0,main="Complete Linkage")
dev.off()
rm(Dmat,com.hclust)

## Use limma to generate moderated t-statistic p-values ######
colnames(a)
contrast <- "sample"
#a <- a_backup
colnames(a)
a <- a[,c(1:3, 10:11)]

treat <- factor(ifelse(grepl(contrast, colnames(a)),"DRUG","CTRL"), levels=c("CTRL","DRUG"))
ctrl <- factor(ifelse(grepl("DMSO",colnames(a)),"DMSO","TREAT"), levels=c("DMSO", "TREAT"))
treat
ctrl

######## LIMMA MODERATED T-STATISTICS #####
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


  # Add gene names to tt for PD1 data
  head(tt)
  
  identifier <- idmapping
  colnames(identifier) <- c("ID", "GeneID")
  identifier <- identifier[,1:2]
  colnames(identifier)
  
  tt <- left_join(x = tt, y = identifier, by = "ID")
  colnames(tt)
  tt <- tt[,c(1,8,2,3,4,5,6,7)]
  colnames(tt)
  rownames(tt) <- tt[,1]
  
  
  # Write hit table
  write.csv(tt, file=paste("./tables/",pre, "_Limma_output_",contrast,"_vs_CTRLs.csv",sep=""), row.names = F)
  
rm(i,sel1,sel2,sel3,sel,b,treat,design,fit,identifier)

  
  
  # Plot & Hits by p.value
  # Define p-value cut-off
  pval <- 0.001
  lab_idx <- log10(pval)+0.4
  # Define minimal logFC
  lfc <- 1.5 # 50% change
  lfc_idx <- log2(lfc) - 0.2
  # find y-axis value
  maxFC <- max(abs(tt$logFC))
  
  pdf(file = paste("plots/",pre,"_hits_",contrast,"_Pvalue.pdf"), width = 12, height = 12)
  sel <- tt$logFC <= -log2(lfc) & tt$P.Value <= pval | tt$logFC >= log2(lfc) & tt$P.Value <= pval
  summary(sel)
  par(mar=c(5,4,4-3,2)+.1, cex=1.8)
  plot(log10(tt$P.Value),tt$logFC,  xlab="log10 P value", ylab=paste(contrast,"_vs_CTRLs -(log fold change)"),
       ylim= c(-(maxFC),maxFC),
       cex=0.8,
       pch=20, col=ifelse(sel, "red", "#22222222"))      
  text(log10(tt[sel,"P.Value"]), tt[sel,"logFC"],adj = c(-0.25,0.55), labels=tt[sel,"GeneID"], col="black", cex=1.0)
  abline(v=log10(pval), lty=2)       
  abline(h=c(-1,1)*log2(lfc), lty=2) 
  legend(pch=20, "topleft", bty="n", legend=paste("Hits: p < ", pval, "; logFC > ",lfc, "-fold"), cex=1.0, col="red")
  text(-3.75, 0.2, labels=paste(lfc, "-fold upregulated  "), cex = 1)
  text(-3.72,-0.2, labels=paste(lfc, "-fold downregulated"), cex = 1)
  text(-0.6,-0.8, labels = paste("P = ",pval), cex = 1)
  dev.off()
 
  
  
  
  # Plot & Hits by adjusted p.value
  # # Define p-value cut-off
  pval <- 0.01
  lab_idx <- log10(pval)+0.35
  # Define minimal logFC
  lfc <- 1.5 # 50% change
  lfc_idx <- log2(lfc) - 0.25
  # find y-axis value
  maxFC <- max(abs(tt$logFC))
  #maxFC <- 2

  pdf(file = paste("plots/",pre,"_hits_",contrast,"_adjPvalue.pdf"), width = 12, height = 12)
  sel <- tt$logFC <= -log2(lfc) & tt$adj.P.Val <= pval | tt$logFC >= log2(lfc) & tt$adj.P.Val <= pval
  summary(sel)
  #plot
  par(mar=c(5,4,4-3,2)+.1, cex=1.8)
  plot(log10(tt$adj.P.Val),tt$logFC,  xlab="log10 adjusted P value", ylab=paste(contrast," - CTRL (log fold change)"),
       ylim= c(-(maxFC),maxFC),
       pch=20, col=ifelse(sel, "red", "#22222222"))
  text(log10(tt[sel,"adj.P.Val"]), tt[sel,"logFC"],adj = c(-0.25,0.55), labels=tt[sel,"GeneID"], col="black", cex=1.0)
  abline(v=log10(pval), lty=2)
  abline(h=c(-1,1)*log2(lfc), lty=2)
  legend(pch=20, "topleft", bty="n", legend=paste("Hits: adj p < ", pval, "; logFC > ",lfc, "-fold"), cex=1.0, col="red")
  text(-2.8, 0.2, labels=paste(lfc, "-fold upregulated  "), cex = 1)
  text(-2.79,-0.2, labels=paste(lfc, "-fold downregulated"), cex = 1)
  text(-0.2,-1, labels = paste("P = ",pval), cex = 1)
  dev.off()

