# Run WGCNA on mCNV data 

#Don't run this script on a macOS as the calculation will be incorcet. 


options(stringsAsFactors=FALSE)
packages <- c("impute","dynamicTreeCut","flashClust","Hmisc",
              "WGCNA","gridExtra","grid","gtable","dplyr")
sapply(packages, require, character.only = TRUE)


#run  WGCNA

outputFolder <- "output/wgcna/"

load("data/GeneAnnotation.rda")

conditions <- c("cortex", "hippocampus")


for (condition in conditions) {
  
  print(condition)
  load(file = paste0(outputFolder, condition, "_star_cqn_noramlize_regress.rdata"))
  datExpr <- t(get(paste0("datExpr_",condition)))
  datMeta <- get(paste0("datMeta_",condition))
  if (!exists("geneAnno")) {
    geneAnno <-  geneAnnoRaw %>% 
      group_by(ensembl_gene_id) %>%
      filter(chromosome_name %in% c(1:22,"X","Y","MT")) %>%   ##  only keep real chromosomes
      filter(transcript_length == max(transcript_length)) %>% ## choose the longest transcript for each gene
      as.data.frame()
    geneAnno <- geneAnno[match(rownames(get(paste0("datExpr_",condition))),geneAnno[,1]),]
    geneAnno <- geneAnno[!is.na(geneAnno$ensembl_gene_id),]
  }
  
  # create subset of datMeta and format it for WGCNA
  # make genotype the first cofactor
  cofactors <- c("Genotype","CNV","Timepoint","Gender","Sequencing.study")
  datMeta1 <- datMeta[,cofactors]
  datMeta1[,] <-  lapply( datMeta1[,] , as.factor)
  datMeta1[,"Genotype"] <- factor(datMeta1[,"Genotype"],levels = c("wt","hemi"))
  
  datMeta1[,] <-  lapply( datMeta1[,] , as.numeric) 
  datMeta1[,"Genotype"]<-datMeta1[,"Genotype"]-1 #0=NEG, 1=POS
  datMeta1 <- datMeta1[,sapply( datMeta1[,] , function(x) length(unique(x))>1) ] #drop cofactos with all the same entries
  
  
   # test that genes don't have too many missing values
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  print(gsg$allOK)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 <- datExpr
    datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  
  
  #correlate genes with traits
  geneSigs=matrix(NA,nrow=ncol(datMeta1),ncol=ncol(datExpr)) # create a vector to hold the data
  rownames(geneSigs)=colnames(datMeta1)
  
  for(i in 1:ncol(geneSigs)) {
    exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
    #Genotype or Dx
    geneSigs[1,i]=cor(exprvec, datMeta1[,1],use="pairwise.complete.obs",method="spearman") 
    for (j in (2:ncol(datMeta1))){
      geneSigs[j,i]=sqrt(max(summary(lm(exprvec~as.factor(datMeta1[,j])))$adj.r.squared,0)) # calculate adjusted R^2s square-root for categorical variables
    }
    if (i%%1000==0){cat('Done for gene...',i,'\n')}
  }
  geneSigs[1,] = numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 
  for (j in (2:ncol(datMeta1))){
    geneSigs[j,] = numbers2colors(as.numeric(geneSigs[j,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  }
  
  
  # Choose a soft-threshold
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  # Plot the results:
  
  par(mfrow = c(1,2))
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.8,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

   
  #set soft threshold
  if (condition == "cortex"){
    softPower = 4
  } else if (condition == "hippocampus") {
    softPower = 5
  } 
  
  #create dissimilarity TOM
  adjacency = adjacency(datExpr, power = softPower, type = "signed",corFnc = "bicor",corOptions = "use='pairwise.complete.obs'")
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1 - TOM
  save(dissTOM, file = paste0(outputFolder,condition,"_dissTOM.Rdata" ))

  # set parameters
  if (condition == "cortex") {
    ds = 2
    mms = 100
    dthresh = 0.1
    
  } else if (condition == "hippocampus") {
    ds = 4
    mms = 160
    dthresh = 0.1
    
  } 
  # create Modules
  createModules(outputFolder,condition)
  
  # Call the hierarchical clustering function
  geneTree = flashClust(as.dist(dissTOM), method = "average");  
  
  tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                      minClusterSize = mms, cutHeight = 0.9999,
                      deepSplit = ds, distM = as.matrix(dissTOM))
  merge <- mergeCloseModules(exprData = datExpr,colors = tree$labels, cutHeight = dthresh)
  mColorh <- cbind(labels2colors(merge$colors),t(geneSigs))
  mLabelh <- c("Merged Colors",rownames(geneSigs))
  
  pdf(paste0(outputFolder,condition,"ModuleDendro.pdf"),height = 10, width = 16) 
  plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = ",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh))
  dev.off()
  
  mergedColors = labels2colors(merge$colors);
  
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  MEList=moduleEigengenes(datExpr, colors = mergedColors,softPower= softPower, nPC=1)
  MEs=MEList$eigengenes
  rownames(MEs) <- rownames(datExpr)
  MEs=orderMEs(MEs)
  moduleColors = mergedColors
  table(moduleColors)
  
  KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")
  
  geneAnnoSubset <- geneAnno[match(rownames(t(datExpr)),geneAnno$ensembl_gene_id),]
  
  geneInfo=as.data.frame(cbind(geneAnnoSubset$ensembl_gene_id,geneAnnoSubset$mgi_symbol,moduleColors, KMEs))  
  colnames(geneInfo)[1]= "Ensembl.Gene.ID"
  colnames(geneInfo)[2]= "GeneSymbol"
  colnames(geneInfo)[3]= "Module"
  write.csv(geneInfo,file=paste0(outputFolder,condition,"geneModules.csv"))
  #save data for later use
  assign(paste0("geneInfo",condition),geneInfo)
  assign(paste0("tree",condition),tree)
  assign(paste0("mColorh",condition),mColorh)
  assign(paste0("mLabelh",condition),mLabelh)
  assign(paste0("merge",condition),merge)
  assign(paste0("KMEs",condition),KMEs)
  assign(paste0("moduleColors",condition),moduleColors)
  assign(paste0("MEList",condition),MEList)
  assign(paste0("MEs",condition),MEs)
  assign(paste0("geneTree",condition),geneTree)
  assign(paste0("datExpr",condition),datExpr)
  rm(list=ls()[grep("dissTOM",ls())])
  save(list=ls()[grep(paste0(condition,"$"),ls())], file=paste0(outputFolder,condition,"_wgcna.Rdata" ))
  
  
}
