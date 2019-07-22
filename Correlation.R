### Co-expression analysis ####
library(corrplot)
library(Hmisc)
library(reshape2)

#The following lines load the data and process it. We use only 113 genes to speed up
targenes=read.delim("Genes_selected.txt", sep = "\t", row.names = 1)
ens2gene=read.delim("Ensembl2gene.tsv", row.names = 2, stringsAsFactors = F)
FPKMdata=read.csv("FPKM_selected.txt", sep = "\t", row.names = 1)
FPKMdata=FPKMdata[row.names(targenes),]
row.names(FPKMdata)=ens2gene[row.names(FPKMdata),]
FPKMdata=t(as.matrix(FPKMdata[rowMeans(FPKMdata)>5,])) # Excludes lowly expressed

## Build the correlation matrix
cor.matrix=rcorr(as.matrix(FPKMdata), type="spearman")
cor.matrix.R=cor.matrix$r
cor.matrix.P=cor.matrix$P

### !! Important! At this step, you should correct for multiple hypothesis testing. For simplicity, we will only use P, but you should always correct P
#Plot
corrplot(cor.matrix.R, p.mat = cor.matrix.P, 
         order = "FPC", #columns and rows sorted by hierarchical clustering
         insig = "blank", sig.level = 0.001,
         addgrid.col=NA,method='square', tl.cex=0.3
) #shows correlations with P>0.001 as blanks

#The following code corrects for multiple hypothesis testing and outputs
uppTriagR=cor.matrix.R; uppTriagR[lower.tri(uppTriagR, diag = T)]=NA;  
uppTriagR=melt(uppTriagR); colnames(uppTriagR)=c("Gene1","Gene2","R"); 
uppTriagP=cor.matrix.P; uppTriagP[lower.tri(uppTriagP, diag = T)]=NA; uppTriagP=melt(uppTriagP); colnames(uppTriagP)=c("Gene1","Gene2","P")
combinedCorr=cbind(uppTriagR,uppTriagP[,"P",drop=F]); combinedCorr=combinedCorr[!is.na(combinedCorr$R),]
combinedCorr[,"FDR"]=p.adjust(combinedCorr$P, method = "fdr") #FDR
combinedCorr=combinedCorr[combinedCorr$FDR<0.05,]
# temp=combinedCorr[,c("Gene2","Gene1","R","P","FDR")]; colnames(temp)= c("Gene1","Gene2","R","P","FDR")
# ToCytoscape=rbind(combinedCorr,temp)
ToCytoscape=data.frame(combinedCorr, stringsAsFactors = F)
row.names(ToCytoscape)=NULL
write.table(ToCytoscape,"Coexpression_113genes.txt", sep = "\t", row.names = F, quote = FALSE)