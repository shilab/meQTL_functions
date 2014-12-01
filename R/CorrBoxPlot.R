#' @export
CorrBoxPlot <-
  function (mEQTL,threshold,expr,genot,visual=FALSE,pdf_file="",crlt=0,cis=TRUE){
    # Inputs:
    #   mEQTL     - Matrix EQTL object with the eQTLs already collected
    #   threshold - FDR cutoff, only those eQTLs with equal or lower threshold will be taken into account
    #   expr      - Transcript expression dataset
    #   genot     - Genotyping dataset, either phased or unphased variants
    #   cis       - If TRUE only cis eQTLs are considered, otherwise trans eQTLS (default TRUE)
    #   visual    - If TRUE the script will display a box plot figure for each eQTL above the threshold
    #
    # Output: A vector with Pearson correlation scores for each eQTL that surpasses the given threshold
    #
    # NOTES:
    # Obviously the original files from which mEQTL object was computed must match on transcript, variants
    # and samples IDs included in expr and genot 
    #
    # expr and genot datafiles are in the matrixEQTL format and can be loaded as:
    # expr = read.table(file_name, header = TRUE, stringsAsFactors = FALSE);
    #
    # R. Armananzas - Last update 12/02/13
    #
    
    
    corr  <- NULL; phenotype <- NULL; genotype <- NULL; phenotype_new<-NULL;
    
    if (cis==TRUE){
      index <- which(mEQTL$cis$eqtls$FDR<=threshold)
      eqtls <- mEQTL$cis$eqtls[index,]
    }
    else{
      index <- which(mEQTL$trans$eqtls$FDR<=threshold)
      eqtls <- mEQTL$trans$eqtls[index,]    
    }
    
    for (i in 1:nrow(eqtls)){
      phenotype[[i]] <- expr[which(expr[,1]==as.character(eqtls$gene[i])),2:ncol(expr)]
      genotype[[i]]  <- genot[which(genot[,1]==as.character(eqtls$snps[i])),2:ncol(genot)]
      corr[i]   <- cor(as.numeric(phenotype[[i]]),as.numeric(genotype[[i]]), use="pairwise.complete.obs")
    }
    
    if (visual){ #Perform the plots
      #There can be three values (unphased) or four (phased)
      if (pdf_file!="")
      {
        pdf(paste('./',pdf_file,sep=""))
        par(mfcol = c(2, 2))
      }
      genotypes <- range(genot[,2:ncol(genot)])[1]:range(genot[,2:ncol(genot)])[2]
      for (i in 1:nrow(eqtls)){
        #Prepare the matrix
        pheno <- as.numeric(phenotype[[i]])
        values <- NULL;
        for (j in 1:length(genotypes)){
          values[[j]] <- pheno[which(genotype[[i]]==genotypes[j])]
        }
        #Plot the boxplots
        #if (length(genotypes)==3){cats=c(0,1,2)}
        #else {cats=c(0,1,2,3,4)}
        cats=seq(0,length(genotypes)-1)
        if (abs(corr[i])>=crlt)
        {
          boxplot(values,boxwex=0.5,ylab=paste(as.character(eqtls$gene[i])," expression"), names=cats,
                  xlab=paste(as.character(eqtls$snps[i])," genotype","\nCorrelation: ",format(corr[i],2),
                             "P-value: ",format(eqtls$pvalue[i],2)," FDR: ",format(eqtls$FDR[i],2)),
                  main=paste(as.character(eqtls$snps[i])," - ",as.character(eqtls$gene[i])))
        }
      }
      if (pdf_file!="")
      {
        dev.off()
      }
    }
    return(corr)
  }