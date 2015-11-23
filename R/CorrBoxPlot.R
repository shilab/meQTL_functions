#' @export
CorrBoxPlot <- function (mEQTL,threshold,expr,genot,visual=FALSE,pdf_file="",crlt=0,cis=TRUE)
{
	# Inputs:
	#   mEQTL     - Matrix EQTL object with the eQTLs already collected
	#   threshold - FDR cutoff, only those eQTLs with equal or lower threshold will be taken into account
	#   expr      - Transcript expression dataset
	#   genot     - Genotyping dataset, either phased or unphased variants
	#   visual    - If TRUE the script will display a box plot figure for each eQTL above the threshold
	#   pdf_file  - Filename for plots
	#   crlt      - Correlation cutoff (default 0)
	#   cis       - If TRUE only cis eQTLs are considered, otherwise trans eQTLS (default TRUE)
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
	# R. Armananzas and Andrew Quitadamo


	corr  <- NULL; phenotype <- NULL; genotype <- NULL;
	index <- getIndex(cis, mEQTL, threshold)
	eqtls <- getEQTLS(cis, mEQTL, index)

	for (i in 1:nrow(eqtls))
	{
		phenotype[[i]] <- expr[which(rownames(expr)==as.character(eqtls$gene[i])),2:ncol(expr)]
		genotype[[i]] <- genot[which(rownames(genot)==as.character(eqtls$snps[i])),2:ncol(genot)]
		corr[i] <- getCorr(phenotype[[i]], genotype[[i]])
	}

	if (visual)
	{
		if (pdf_file!="")
		{
			pdf(paste('./',pdf_file,sep=""))
			par(mfcol = c(2, 2))
		}
		genotypes <- range(genot[,2:ncol(genot)])[1]:range(genot[,2:ncol(genot)])[2]
		for (i in 1:nrow(eqtls))
		{
			#Prepare the matrix
			pheno <- as.numeric(phenotype[[i]])
			values <- list();
			for (j in 1:length(genotypes))
			{
				values[[j]] <- pheno[which(genotype[[i]]==genotypes[j])]
			}
			#Plot the boxplots
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

getIndex <- function(cis, mEQTL, threshold)
{
	if (cis==TRUE)
   	{
       return(which(mEQTL$cis$eqtls$FDR<=threshold))
	}
	else
	{
		return(which(mEQTL$trans$eqtls$FDR<=threshold))
	}
}

getEQTLS <- function(cis, mEQTL, index)
{
	if (cis==TRUE)
	{
		return(mEQTL$cis$eqtls[index,])
	}
	else
	{
		return(mEQTL$trans$eqtls[index,])		
	}
}

getCorr <- function(phenotype, genotype)
{
	return(cor(as.numeric(phenotype),as.numeric(genotype), use="pairwise.complete.obs"))
}
