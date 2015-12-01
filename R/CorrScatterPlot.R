#' @export
CorrScatterPlot <- function (mEQTL,threshold,expr,genot,visual=TRUE,pdf_file="",crlt=0,cis=TRUE)
{
  # Inputs:
  #   mEQTL     - Matrix EQTL object with the eQTLs already collected
  #   threshold - FDR cutoff, only those eQTLs with equal or lower threshold will be taken into account
  #   expr      - Expression filename
  #   genot     - Genotype filename
  #   cis       - If TRUE only cis eQTLs are considered, otherwise trans eQTLS (default TRUE)
  #   visual    - If TRUE the script will display a box plot figure for each eQTL above the threshold
  #
  # Output: A vector with Pearson correlation scores for each eQTL that surpasses the given threshold
  #
  # NOTES:
  # Obviously the original files from which mEQTL object was computed must match on transcript, variants
  # and samples IDs included in expr and genot 
  #
  # R. Armananzas & A. Quitadamo
  
	expr <- getFileData(expr)
	genot <- getFileData(genot)

	index <- getIndex(cis, mEQTL, threshold)
	eqtls <- getEQTLS(cis, mEQTL, index)
  
	phenotype <- getEQTLPhenotypes(eqtls, expr)
	genotype <- getEQTLGenotypes(eqtls, genot)

	corr <- mapply(getCorr, phenotype, genotype)

	if (visual)
	{ #Perform the plots
		if (pdf_file!="")
		{
			pdf(pdf_file)
			par(mfcol = c(2, 2))
		}

		for (i in 1:nrow(eqtls))
		{
			#Prepare the matrix
			geno <- as.numeric(genotype[[i]])
			pheno <- as.numeric(phenotype[[i]])
			if (abs(corr[i])>=crlt)
			{
				plot(pheno, geno, xlab="miRNA Expression", ylab="mRNA Expression", main = paste(as.character(eqtls$snps[i])," genotype","\nCorrelation: ",format(corr[i],2),"P-value: ",format(eqtls$pvalue[i],2)," FDR: ",format(eqtls$FDR[i],2)))
			}
		}
	}
	if (pdf_file!="")
	{
		dev.off()
	}
	return(corr)
}
