setModel <-function(model)
{
    if (model=="linear")
    {
       return(modelLINEAR)
    }
    else if (model=="anova")
    {
       return(modelANOVA)
    }
    else if (model=="linear_cross")
    {
       return(modelLINEAR_CROSS)
    }
	else
	{
		stop('Model must be linear, anova or linear_cross');
	}
}

#' @export
mxeqtl <-
function(snp_file,snp_location,expr_file,expr_location,cis_output_file,
         cis_pval,covariates="",trans_output_file="", trans_pval=0, 
         model="linear", MAF=0, cis_dist=1e6, qq="", missing="NA", sep="\t")
{
	# Matrix eQTL function based on the sample code by Andrey A. Shabalin
	# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

	useModel = setModel(model)

	# Genotype file name
	SNP_file_name = snp_file
	snps_location_file_name = snp_location

	# Gene expression file name
	expression_file_name = expr_file
	gene_location_file_name = expr_location

	# Covariates file name
	if (covariates!="")
	{
		covariates_file_name = covariates;
		cat(covariates_file_name,'\n');
	}
	else
	{
		covariates_file_name = character();
	}
	# Output file name
	output_file_name = tempfile();
	output_file_name_cis = cis_output_file
	if (trans_output_file!="")
	{
		output_file_name_tra = trans_output_file;
	}
	else
	{
		output_file_name_tra = "onlyTRANSresults.txt"
	}

	# Set p-value levels
	pvOutputThreshold_cis = cis_pval;
	pvOutputThreshold_tra = trans_pval

	#TODO: Create option for errorCovariance
	errorCovariance = numeric();

	# Distance for local gene-SNP pairs
	cisDist = cis_dist;

	## Load genotype data
	snps = SlicedData$new();
	snps$fileDelimiter = sep;
	snps$fileOmitCharacters = missing; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	snps$LoadFile(SNP_file_name);

	if (MAF>0)
	{
		maf.list = vector('list', length(snps))
		for(sl in 1:length(snps))
		{
			slice = snps[[sl]];
			maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
			maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
		}
		maf = unlist(maf.list)
	}

	cat('SNPs before filtering:',nrow(snps),'\n')

	if (MAF>0)
	{
		snps$RowReorder(maf>MAF);
		cat('SNPs after filtering:',nrow(snps),'\n')
	}

	snps$SaveFile("meQTL_filtered_input")

	# Load gene expression data
	gene = SlicedData$new();
	gene$fileDelimiter = sep;
	gene$fileOmitCharacters = missing; # denote missing values;
	gene$fileSkipRows = 1;          # one row of column labels
	gene$fileSkipColumns = 1;       # one column of row labels
	gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	gene$LoadFile(expression_file_name);

	# Load covariates
	cvrt = SlicedData$new();
	cvrt$fileDelimiter = sep;
	cvrt$fileOmitCharacters = missing; # denote missing values;
	cvrt$fileSkipRows = 1;          # one row of column labels
	cvrt$fileSkipColumns = 1;       # one column of row labels
	if(length(covariates_file_name)>0)
	{
		cvrt$LoadFile(covariates_file_name);
	}

	# Load the genotype and expression positions
	snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
	genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

	#Call MatrixEQTL function
	me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	output_file_name     = output_file_name_tra,
	pvOutputThreshold     = pvOutputThreshold_tra,
	useModel = useModel,
	errorCovariance = errorCovariance,
	verbose = TRUE,
	output_file_name.cis = output_file_name_cis,
	pvOutputThreshold.cis = pvOutputThreshold_cis,
	snpspos = snpspos,
	genepos = genepos,
	cisDist = cisDist,
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = TRUE,
	noFDRsaveMemory = FALSE);

	cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
	cat('Detected ',me$cis$neqtls,' local eQTLs:', '\n');
	cat('Detected ',me$trans$neqtls,' distant eQTLs:', '\n');

	## Plot the Q-Q plot of local and distant p-values
	if (qq!="")
	{
		pdf(qq)
		plot(me)
		dev.off()
	}
	else
	{
		plot(me)
	}

	return(me)
}
