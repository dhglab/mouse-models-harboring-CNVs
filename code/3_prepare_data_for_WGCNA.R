#5/12/16
#Aaron Gordon
#Run data normalization and QC on mCNV data

# set up R and some variables ----------------------------
setwd("~/Aaron/Projects/mCNVs/")
options(stringsAsFactors = FALSE)
library(ggplot2)
library(WGCNA)
library(plyr)
library(cqn)
library(tidyverse)
library(stringr)
library(sva)
rm(list = ls())
outputfile_counter = 1
outputFolder = "output/wgcna"
dir.create(outputFolder,showWarnings = F,recursive = T)

outlierAnalysis <- function(datExpr,datMeta,title,color.by){
  sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
  netsummary <- fundamentalNetworkConcepts(normadj); 
  K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
  C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
  outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
  print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
  outlierIDs <- (rep(" ", length(outliers))); outlierIDs[outliers] <- levels(datMeta$ID2)[as.numeric(datMeta$ID2[outliers])]
  plot(Z.K, col = as.numeric(as.factor(datMeta[,color.by])), pch=19, main=paste("Outlier detection for",title,"\n color = ",color.by),xaxt="n",xlab=NA, ylab="Network connectivity (z score)");text(1:length(Z.K),Z.K,label=outlierIDs,pos=3,cex=0.6);abline(h=-2, lty=2)
  #axis(1, at=1:123, ,labels=datMeta$Day.grouped,las=2)
  return (outliers)
}



#Load Expression, Meta and QC data  and format them ----------------------------------------
load(file = "data/compiledData.rdata") # compiled using 0.2_compile_data.R

datExprRaw <- datExpr
datMetaRaw <- datMeta
datQCRaw <- datQC

datMetaRaw$Sequencing.study <- make.names(datMetaRaw$Sequencing.study)
datMetaRaw$Timepoint[is.na(datMetaRaw$Timepoint)] <- (-1)
datMetaRaw$RIN.value[is.na(datMetaRaw$RIN.value)] <- (-1)
datMetaRaw$CNV <- levels(datMetaRaw$CNV)[datMetaRaw$CNV]
datMetaRaw$CNV <- as.factor(ifelse(datMetaRaw$Sequencing.study == "IBP.ANGF.study" & datMetaRaw$Genotype != "hemi", 
                                   paste(datMetaRaw$CNV, "(homo)"), 
                                   datMetaRaw$CNV))



# subset data --------------------
#only keep hemizygous and matching wt samples

samples_to_keep <- (datMeta$CNV %in% c("1q21","22q11") | grepl("2|4|6",datMeta$Sequencing.study )) & 
  datMeta$Tissue %in% c("cortex", "hippocampus")

datMeta <- subset(datMetaRaw, samples_to_keep)
datMeta   <-  droplevels(datMeta)

conditions = unique(datMeta$Sequencing.study)





#  Filter very low epression genes  --------------------
# genes with less than 10 counts in 60% of  samples 

pres <- apply(datExprRaw>10,1,sum)
idx <- which(pres > 0.8*dim(datExprRaw)[2])
datExprRawFilt <- datExprRaw[idx,]

datMeta_subset <- datMeta


for (condition in conditions) {
 
  #  subset data
  datMeta <- subset(datMeta_subset, Sequencing.study == condition & Genotype != "homo")
  datMeta <- droplevels(datMeta)
  datExpr <- datExprRawFilt[, colnames(datExprRawFilt) %in% rownames(datMeta)]
  datQC <- datQCRaw[rownames(datQCRaw) %in% rownames(datMeta),]
  
  # (4) Get Gene Annotaion --------------------------------------------
  # library(biomaRt)
  # getinfo <- c("ensembl_gene_id","mgi_symbol","chromosome_name","strand","start_position", "end_position","gene_biotype","transcript_length","percentage_gc_content")
  # mouseMart <- useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl") #change if working with differnet spieces
  # geneAnnoRaw <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=rownames(datExprRaw),mart=mouseMart)
  # save(geneAnnoRaw,file="./data/GeneAnnotation.rda")
  
  load("./data/GeneAnnotation.rda")
  geneAnno <-  geneAnnoRaw %>% 
    group_by(ensembl_gene_id) %>%
    filter(chromosome_name %in% c(1:22,"X","Y","MT")) %>%   ##  only keep real chromosomes
    filter(transcript_length == max(transcript_length)) %>% ## choose the longest transcript for each gene
    as.data.frame()
  geneAnno <- geneAnno[match(rownames(datExpr),geneAnno[,1]),]
  geneAnno <- geneAnno[!is.na(geneAnno$ensembl_gene_id),]
  datExpr <- datExpr[match(geneAnno$ensembl_gene_id,rownames(datExpr)),]
  
  
  
  
  
  # (5) Normalize ------------------------------------------------------
  
  
  datExprCQN <- cqn(datExpr,x = geneAnno$percentage_gc_content, lengths = geneAnno$transcript_length,sqn = F ,verbose = T)
  RPKM.cqn <- datExprCQN$y + datExprCQN$offset
  datExpr.htg <- RPKM.cqn
  
  
  # (6) Remove Outliers --------------------------------------------------------------------
 
  
  outliers <- outlierAnalysis(datExpr,datMeta,condition,"CNV")
  datMeta<-datMeta[!outliers,]
  datExpr.htg<-datExpr.htg[,!outliers]
  datQC<-datQC[!outliers,]
  
  
  ## get the PCs of the sequencing statistics ------------------------------------------------------------------
  datQCcolstoRemove <-  grep("CATEGORY|ESTIMATED_LIBRARY_SIZE|pair",colnames(datQC),ignore.case = T)
  datQCforPCA <- datQC[,-c(datQCcolstoRemove)]
  large_qc<- which(apply(datQCforPCA,2,mean)>10^4)
  datQCforPCA[,large_qc] = log10(datQCforPCA[,large_qc])
  PC.datQC <- prcomp(na.omit(t(scale((datQCforPCA),scale=F))), center=T)
  varexp <- (PC.datQC$sdev)^2 / sum(PC.datQC$sdev^2)
  topPC.datQC <- PC.datQC$rotation[,1:5]
  colnames(topPC.datQC) <- c("SeqPC1","SeqPC2" ,"SeqPC3","SeqPC4","SeqPC5")
  
  datMeta <- cbind(datMeta,topPC.datQC)
  

  # Regress out covaritaes ----------------------------------------------------
  datExpr_preRgression <- datExpr

  datMetaforRegress <- datMeta
  datMetaforRegress$Timepoint <- as.character(datMetaforRegress$Timepoint)
  datMetaforRegress$Timepoint[is.na(datMetaforRegress$Timepoint)] <- -1 
  datMetaforRegress$Timepoint <- factor(datMetaforRegress$Timepoint)
  reg_formula <- as.formula("~ Tissue + Genotype + SeqPC1 + SeqPC2 + SeqPC3 + SeqPC4 + SeqPC5")
  if(condition=="Study.2"){
    reg_formula <- as.formula("~ Genotype + SeqPC1 + SeqPC2")
    
  } else if (length(unique(datMeta$Tissue))<2){
    reg_formula <- as.formula("~ Genotype + SeqPC1 + SeqPC2 + SeqPC3 + SeqPC4 + SeqPC5")
    
  }
  
  X =  model.matrix(reg_formula, data = datMetaforRegress) #remove Technical covariates
  Y = datExpr_preRgression
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  
  cols_to_regress <- grep("SeqPC",colnames(X))
  to_regress = (as.matrix(X[,cols_to_regress]) %*% (as.matrix(beta[cols_to_regress,])))  #Timepoint+Gender+CNV+regress out seqPC1-5
  datExpr = datExpr_preRgression - t(to_regress)
  
  
  assign(paste0("datExpr_",condition),datExpr)
  assign(paste0("datMeta_",condition),datMeta)
  assign(paste0("datQC_",condition),datQC)
  save(list=paste(c("datExpr","datMeta","datQC"),condition,sep="_"),file=paste0(file=paste0(outputFolder,condition,"_star_cqn_noramlize_regress.rdata")))
}

# (10) combine all studies and split by Tissue  ------------------------------------------------------
for (conditon in conditions){
  load(file=paste0(file=paste0(file=paste0(outputFolder,condition,"_star_cqn_noramlize_regress.rdata"))))
}
datExprAll <- sapply(conditions, function(x) get(paste0("datExpr_",x)))
datExprAll <- do.call("cbind",datExprAll)
datMetaAll <- lapply(conditions, function(x) get(paste0("datMeta_",x)))
datMetaAll <- do.call("rbind",datMetaAll)
stopifnot(all(rownames(datMetaAll)==colnames(datExprAll)))



# egress out covariates per tissue ---------------------------------
for (condition2 in unique(datMetaAll$Tissue)){
  
  cat(condition2)
  cat("\n")
  
  datMeta <- subset(datMetaAll, datMetaAll$Tissue==condition2)
  
  datExpr <- datExprAll[,colnames(datExprAll) %in% rownames(datMeta)]
  datMetaforRegress <- datMeta
  datMetaforRegress$Timepoint <- as.character(datMetaforRegress$Timepoint)
  datMetaforRegress$Timepoint[is.na(datMetaforRegress$Timepoint)] <- -1 
  datMetaforRegress$Timepoint <- factor(datMetaforRegress$Timepoint)
  datMetaforRegress <- droplevels(datMetaforRegress)
  reg_formula <- as.formula("~ Genotype + CNV ")
  
  X =  model.matrix(reg_formula, data = datMetaforRegress) #remove  covariates
  Y = datExpr
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  
  cols_to_regress <- grep("Intercept|Genotype",colnames(X),invert = T)
  to_regress = (as.matrix(X[,cols_to_regress]) %*% (as.matrix(beta[cols_to_regress,])))  
  datExpr = datExpr - t(to_regress)
  
  # Remove outliers from combined data
  outliers2 <- outlierAnalysis(datExpr,datMeta,condition2, "Sequencing.study")
  
 
  datMeta <- datMeta[!outliers2,]
  datExpr <- datExpr[,!outliers2]  
  
  cat("Saving results\n\n")
  assign(paste0("datExpr_",condition2),datExpr)
  assign(paste0("datMeta_",condition2),datMeta)
  assign(paste0("datQC_",condition2),datQC)
  save(list=paste(c("datExpr","datMeta"),condition2,sep="_"),file=paste0(file=paste0(outputFolder,condition2,"_star_cqn_noramlize_regress.rdata")))
}



