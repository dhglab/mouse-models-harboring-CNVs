#Run limma voom 

#(0) setup and load data ------------------------------------------------------

options(stringsAsFactors = FALSE)
library(limma)
library(tidyverse)
library(reshape)
library(edgeR)
base_folder <- "output/DE/"


load("data/compiledData.rdata")
load("data/GeneAnnotation.rda")


geneAnno <-  geneAnnoRaw %>% 
  filter(chromosome_name %in% c(1:22,"X","Y","MT")) %>%   ##  only keep real chromosomes
  group_by(ensembl_gene_id) %>%
  filter(transcript_length == max(transcript_length)) %>% ## choose the longest transcript for each gene
  as.data.frame()

datMetaRaw <- datMeta
datExprRaw <- datExpr
load(file = paste0(base_folder,"all_norm_data.rdata"))
tissues <- unique(datMetaAll$Tissue)
CNVs <- sort(as.character(unique(datMetaAll$CNV)))


# run limma voom ------------------------------------------------------------
cnv_limmaVoom <- list()
Ns <- data.frame(CNV.Tissue = character(), n = character())
for (tissue in tissues) {
  cnv_limmaVoom[[tissue]] <- list()
  for (cnv in CNVs) {
    datMeta <- datMetaAll %>%
      mutate(row.name = rownames(.)) %>%
      filter(CNV == cnv & Tissue == tissue) %>%
      droplevels() %>%
      mutate(Genotype = relevel(.$Genotype, "wt"))
    
    datExpr <- datExprAll
    tExprRaw <- datExprRaw[rownames(datExpr),datMeta$row.name]
    
    datExpr_DGE <- DGEList(tExprRaw)
    datExpr1 <- calcNormFactors(datExpr_DGE, method = "TMM")
     
    # build model for DE
    covars <- c("Gender", paste0("SeqPC",1:5))
    if (cnv == "15q13" & tissue == "hippocampus") {
      covars <- c("Gender")
      
    }
    design1 <- "~ Genotype"
    for (covar in covars) {
      if (length(unique(datMeta[,covar])) >= 2) {
        design1 <- paste(design1,covar,sep = " + ")
      }
    }
    
    cat(paste("The formula used for DE for", cnv, tissue, "is:" ,design1,"\n"))
    mod1 = model.matrix(as.formula(design1),data = datMeta)
    colnames(mod1) = make.names(colnames(mod1))
    
    voomExpr <- voom(datExpr1,mod1) 
    fit <- lmFit(voomExpr,mod1)
    efit <- eBayes(fit)
    res <- topTable(efit, coef=2 ,number = Inf, confint = T) %>% 
      mutate(ensembl_gene_id = rownames(.)) %>%
      left_join(geneAnno) # anotate genes 
    cnv_limmaVoom[[tissue]][[cnv]] <- res      
  }
}

save(cnv_limmaVoom,Ns,file=file.path(base_folder,"voom_results.rdata"))


