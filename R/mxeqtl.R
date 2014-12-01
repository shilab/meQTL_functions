#' @export
mxeqtl <-
function(snp_file,snp_location,expr_file,expr_location,cis_output_file,
         cis_pval,covariates="",trans_output_file="", trans_pval=0, 
         model="linear", MAF=0, cis_dist=1e6, missing="NA")
{
# Matrix eQTL function based on the sample code by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
# modelANOVA, modelLINEAR, or modelLINEAR_CROSS
if (model=="linear")
{
	useModel = modelLINEAR 
}
else if (model=="anova")
{
	useModel = modelANOVA
}
else if (model=="linear_cross")
{
	useModel = modelLINEAR_CROSS
}

missing_data=missing;

# Genotype file name
SNP_file_name = snp_file;
snps_location_file_name = snp_location;

# Gene expression file name
expression_file_name = expr_file;
gene_location_file_name = expr_location;


# Covariates file name
# Set to character() for no covariates
if (covariates!="")
{
	#covariates_file_name = paste(base.dir,covaraites,sep="");
  covariates_file_name = covariates;
}
else
{
	covariates_file_name = character();
}
# Output file name
output_file_name = tempfile();
output_file_name_cis = cis_output_file;
if (trans_output_file!="")
{
  output_file_name_tra = trans_output_file;
}
else
{
	output_file_name_tra = "./onlyTRANSresults.txt"
}

# Only associations significant at this level will be saved
# Only cis-associations are computed
pvOutputThreshold_cis = cis_pval;
pvOutputThreshold_tra = trans_pval

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = cis_dist;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = missing; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

if (MAF>0)
{
		maf.list = vector('list', length(snps))
		for(sl in 1:length(snps)) {
		  slice = snps[[sl]];
		  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
		  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
		}
		maf = unlist(maf.list)
}

cat('SNPs before filtering:',nrow(snps),'\n')
# snps$RowReorderSimple(maf>0.05);
if (MAF>0)
{
	snps$RowReorder(maf>MAF);
	cat('SNPs after filtering:',nrow(snps),'\n')
}

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = missing; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = missing; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

## Remember this: 
## Set pvOutputThreshold = 0 and pvOutputThreshold.cis > 0 to perform eQTL analysis for local gene-SNP pairs only. 
## Local associations significant at pvOutputThreshold.cis level will be recorded in output_file_name.cis.

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
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
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected ',me$cis$neqtls,' local eQTLs:', '\n');
#show(me$cis$eqtls)
cat('Detected ',me$trans$neqtls,' distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values
plot(me)

## Auxiliar plot + change legends
#plot(GD462_CEU_OldGen,xmin=1e-09,ymin=1e-165,main=NULL)
#legend("topleft",legend=c("Cis p-values","Trans p-values","diagonal"),
#       text.col=c("red","blue","grey"),col=c("red","blue","grey"),
#       pch=c(19,19,NA),lty=c(1,1,1),cex=1.06)
return(me)
}
