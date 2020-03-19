# remove outliers and calculate sequencing PCs

# set up R and some variables ----------------------------

options(stringsAsFactors = FALSE)
library(data.table)
library(ggplot2)
library(WGCNA)
library(plyr)
library(cqn)
library(tidyverse)
library(stringr)
library(gridExtra)
library(sva)

###  ### ### ### ### ### ### ### ### ### ### ### ### ### 
#change the working directory to the GitHub directory
setwd("~/repo/mouse-models-harboring-CNVs/")
###  ### ### ### ### ### ### ### ### ###### ### ### ###



outputFolder ="output/DE/"
dir.create(outputFolder, showWarnings = F, recursive = T)

outlierAnalysis <- function(datExpr,datMeta,title,color.by){
  sdout <- 2; 
  normadj <- (0.5 + 0.5 * bicor(datExpr, use = 'pairwise.complete.obs'))^2
  netsummary <- fundamentalNetworkConcepts(normadj); 
  K <- netsummary$Connectivity
  Z.K <- (K - mean(K)) / sqrt(var(K))
  C <- netsummary$ClusterCoef
  Z.C = (C - mean(C)) / sqrt(var(C))
  outliers <- (Z.K > mean(Z.K) + sdout * sd(Z.K)) | (Z.K < mean(Z.K) - sdout * sd(Z.K))
  print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep = "")); 
  print(colnames(datExpr)[outliers]); print(table(outliers))
  outlierIDs <- (rep(" ", length(outliers))); outlierIDs[outliers] <- levels(datMeta$ID2)[as.numeric(datMeta$ID2[outliers])]
  plot(Z.K, col = as.numeric(as.factor(datMeta[,color.by])), pch = 19, 
       main = paste("Outlier detection for",title,"\n color = ",color.by),
       xaxt = "n",xlab = NA, ylab = "Network connectivity (z score)");
  text(1:length(Z.K), Z.K, label = outlierIDs, pos = 3, cex = 0.6);
  abline(h = -2, lty = 2)
  #axis(1, at=1:123, ,labels=datMeta$Day.grouped,las=2)
  return(outliers)
}

#Load Expression, Meta and QC data and format them ----------------------------------------
load(file = "data/compiledData.rdata") 
rownames(datMeta) <- datMeta$matchID
datMetaRaw <- datMeta
datExprRaw <- datExpr[,match(rownames(datMetaRaw), colnames(datExpr))]

datQCRaw <- datQC[match(rownames(datMetaRaw),rownames(datQC)),]

datMetaRaw$Sequencing.study <- make.names(datMetaRaw$Sequencing.study)
datMetaRaw$Timepoint[is.na(datMetaRaw$Timepoint)] <- (-1)
datMetaRaw$RIN.value[is.na(datMetaRaw$RIN.value)] <- (-1)
datMetaRaw$CNV <- levels(datMetaRaw$CNV)[datMetaRaw$CNV]
datMetaRaw$CNV <- as.factor(ifelse(datMetaRaw$Sequencing.study == "IBP.ANGF.study" & datMetaRaw$Genotype != "hemi", 
                                   paste(datMetaRaw$CNV, "(homo)"), 
                                   datMetaRaw$CNV))

stopifnot((colnames(datExprRaw) == rownames(datMetaRaw)))
stopifnot(all(rownames(datQCRaw) == rownames(datMetaRaw)))

# Filter very low epression genes  --------------------
# genes with less than 10 counts in 80% of  samples
pres <- apply(datExprRaw>10,1,sum)
idx <- which(pres > 0.8*dim(datExprRaw)[2])
datExprRawFilt <- datExprRaw[idx,]


# remove outliers and calculate sequencing PCs per study ----------------------------------------------------------
conditions = unique(datMetaRaw$Sequencing.study)
for (condition in conditions) {
  cat(paste0("\n",condition,"\n"))
  outputfile_counter <- 2
  #    subset data
  datMeta <- subset(datMetaRaw, Sequencing.study == condition)
  datMeta <- droplevels(datMeta)
  datExpr <- datExprRawFilt[, colnames(datExprRawFilt) %in% rownames(datMeta)]
  datQC <- datQCRaw[rownames(datQCRaw) %in% rownames(datMeta),]
  
  #  Get Gene Annotaion --------------------------------------------
  if (!file.exists("data/GeneAnnotation.rda")) {
    library(biomaRt)
    getinfo <- c("ensembl_gene_id","mgi_symbol","chromosome_name","strand","start_position", "end_position","gene_biotype","transcript_length","percentage_gc_content")
    mouseMart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") 
    geneAnnoRaw <- getBM(attributes = getinfo,filters = c("ensembl_gene_id"),values = rownames(datExprRaw),mart = mouseMart)
    save(geneAnnoRaw,file = "data/GeneAnnotation.rda")
  } else {
    load("data/GeneAnnotation.rda")
  }
  geneAnno <-  geneAnnoRaw %>% 
    group_by(ensembl_gene_id) %>%
    filter(chromosome_name %in% c(1:22,"X","Y","MT")) %>%   ##  only keep real chromosomes
    filter(transcript_length == max(transcript_length)) %>% ## choose the longest transcript for each gene
    as.data.frame()
  geneAnno <- geneAnno[match(rownames(datExpr),geneAnno[,1]),]
  geneAnno <- geneAnno[!is.na(geneAnno$ensembl_gene_id),]
  datExpr <- datExpr[match(geneAnno$ensembl_gene_id,rownames(datExpr)),]
  
  # CQN Normalization ------------------------------------------------------
  datExpr_preNorm <- datExpr
  datExprCQN <- cqn(datExpr,x=geneAnno$percentage_gc_content, lengths= geneAnno$transcript_length,sqn=F ,verbose = T)
  RPKM.cqn <- datExprCQN$y + datExprCQN$offset
  datExpr.htg<-RPKM.cqn
  
  #  Remove Outliers --------------------------------------------------------------------
  outliers <- outlierAnalysis(datExpr, datMeta, "all samples", "CNV")
  datMeta <- datMeta[!outliers, ]
  datExpr.htg <- datExpr.htg[, !outliers]
  datQC <- datQC[!outliers,]
  
  ## get the PCs of the sequencing statistics ---------------------------------------------------------------------------------
  datQCcolstoRemove <-  grep("CATEGORY|ESTIMATED_LIBRARY_SIZE|pair",colnames(datQC),ignore.case = T)
  datQCforPCA <- datQC[,-c(datQCcolstoRemove)]
  large_qc<- which(apply(datQCforPCA,2,mean)>10^4)
  datQCforPCA[,large_qc] = log10(datQCforPCA[,large_qc])
  PC.datQC <- prcomp(na.omit(t(scale((datQCforPCA),scale=F))), center=T)
  varexp <- (PC.datQC$sdev)^2 / sum(PC.datQC$sdev^2)
  topPC.datQC <- PC.datQC$rotation[,1:5]
  colnames(topPC.datQC) <- c("SeqPC1","SeqPC2" ,"SeqPC3","SeqPC4","SeqPC5")
  
  datMeta <- cbind(datMeta,topPC.datQC) # add sequencing PCs to meta data
  
  assign(paste0("datExpr_",condition),datExpr.htg)
  assign(paste0("datMeta_",condition),datMeta)
  assign(paste0("datQC_",condition),datQC)
  save(list=paste(c("datExpr","datMeta","datQC"),condition,sep="_"),file=file.path(outputFolder,paste0(condition,"_star_cqn_noramlize_regress.rdata")))
}
# 
# combine all studies  ------------------------------------------------------
for (condition in conditions){
  load(file=paste0(outputFolder,condition,"_star_cqn_noramlize_regress.rdata"))
}
datExprAll <- sapply(conditions, function(x) get(paste0("datExpr_",x)))
datExprAll <- do.call("cbind",datExprAll)

datMetaAll <- lapply(conditions, function(x) get(paste0("datMeta_",x)))
datMetaAll <- do.call("rbind",datMetaAll)
stopifnot(all(rownames(datMetaAll)==colnames(datExprAll)))
save(datExprAll, datMetaAll, file = file.path(outputFolder,"all_norm_data.rdata"))
