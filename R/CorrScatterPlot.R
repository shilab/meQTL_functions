#' @export
CorrScatterPlot <- function (mEQTL,threshold,expr,genot,visual=TRUE,cis=TRUE)
{
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
  # R. Armananzas & A. Quitadamo
  
    index <- getIndex(cis, mEQTL, threshold)
	eqtls <- getEQTLS(cis, mEQTL, index)
  
    phenotype <- getEQTLPhenotypes(eqtls, expr)
	genotype <- getEQTLGenotypes(eqtls, genot)

    corr <- mapply(getCorr, phenotype, genotype)

	if (visual)
	{ #Perform the plots
    	for (i in 1:nrow(eqtls))
		{
      	#Prepare the matrix
			geno <- as.numeric(genotype[[i]])
			pheno <- as.numeric(phenotype[[i]])
			plot(pheno, geno, xlab="miRNA Expression", ylab="mRNA Expression", main = paste(as.character(eqtls$snps[i])," genotype","\nCorrelation: ",format(corr[i],2),"P-value: ",format(eqtls$pvalue[i],2)," FDR: ",format(eqtls$FDR[i],2)))
		}
	}
	return(corr)
}
