\name{mxeqtl}
\alias{mxeqtl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix eQTL function
%%  ~~function to do ... ~~
}
\description{Functionalized version of Matrix eQTL
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mxeqtl(snp_file, snp_location, expr_file, expr_location, cis_output_file, cis_pval, covariates = "", trans_output_file = "", trans_pval = 0, model = "linear", MAF = 0, cis_dist = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp_file}{
    File that contains the genotype data. Required.
}
  \item{snp_location}{
    File that contains genotype genetic locations. Required.
}
  \item{expr_file}{
File that contains the expression data. Required.
}
  \item{expr_location}{
File that contains expression location data. Required.
}
  \item{cis_output_file}{
    Name of file cis eQTLs will be output to. Required.
}
  \item{cis_pval}{
    P-value cutoff for cis eQTLs. Required.
}
  \item{covariates}{
    File containing covariate data. Optional.
}
  \item{trans_output_file}{
    Name of file trans eQTLs will be output to. Optional. 
}
  \item{trans_pval}{
    P-value cutoff for trans eQTLs. Optional.
}
  \item{model}{
    Model used for analysis. Default is modelLINEAR.
    Can be set to modelANOVA and modelLINEAR_CROSS.
    Optional.
}
  \item{MAF}{
    MAF value for filtering. Optional. 
}
  \item{cis_dist}{
    Maximum distance for cis eQTLs. Default is 1 megabase (1e6 pb).
}
}
\details{
    Optional arguments need to be called with their name.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrew Quitadamo
}
\note{
Based off the example code of Andrey A. Shabalin.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (snp_file, snp_location, expr_file, expr_location, cis_output_file, 
    cis_pval, covariates = "", trans_output_file = "", trans_pval = 0, 
    model = "linear", MAF = 0, cis_dist = 1e+06) 
{
    base.dir <- "./"
    if (model == "linear") {
        useModel = modelLINEAR
    }
    else if (model == "anova") {
        useModel = modelANOVA
    }
    else if (model == "linear_cross") {
        useModel = modelLINEAR_CROSS
    }
    SNP_file_name = paste(base.dir, snp_file, sep = "")
    snps_location_file_name = paste(base.dir, snp_location, sep = "")
    expression_file_name = paste(base.dir, expr_file, sep = "")
    gene_location_file_name = paste(base.dir, expr_location, 
        sep = "")
    if (covariates != "") {
        covariates_file_name = paste(base.dir, covaraites, sep = "")
    }
    else {
        covariates_file_name = character()
    }
    output_file_name = tempfile()
    output_file_name_cis = paste(base.dir, cis_output_file, sep = "")
    if (trans_output_file != "") {
        output_file_name_tra = paste(base.dir, trans_output_file, 
            sep = "")
    }
    else {
        output_file_name_tra = "./onlyTRANSresults.txt"
    }
    pvOutputThreshold_cis = cis_pval
    pvOutputThreshold_tra = trans_pval
    errorCovariance = numeric()
    cisDist = cis_dist
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile(SNP_file_name)
    if (MAF > 0) {
        maf.list = vector("list", length(snps))
        for (sl in 1:length(snps)) {
            slice = snps[[sl]]
            maf.list[[sl]] = rowMeans(slice, na.rm = TRUE)/2
            maf.list[[sl]] = pmin(maf.list[[sl]], 1 - maf.list[[sl]])
        }
        maf = unlist(maf.list)
    }
    cat("SNPs before filtering:", nrow(snps), "\n")
    if (MAF > 0) {
        snps$RowReorder(maf > MAF)
        cat("SNPs after filtering:", nrow(snps), "\n")
    }
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile(expression_file_name)
    cvrt = SlicedData$new()
    cvrt$fileDelimiter = "\t"
    cvrt$fileOmitCharacters = "NA"
    cvrt$fileSkipRows = 1
    cvrt$fileSkipColumns = 1
    if (length(covariates_file_name) > 0) {
        cvrt$LoadFile(covariates_file_name)
    }
    snpspos = read.table(snps_location_file_name, header = TRUE, 
        stringsAsFactors = FALSE)
    genepos = read.table(gene_location_file_name, header = TRUE, 
        stringsAsFactors = FALSE)
    me = Matrix_eQTL_main(snps = snps, gene = gene, output_file_name = output_file_name_tra, 
        pvOutputThreshold = pvOutputThreshold_tra, useModel = useModel, 
        errorCovariance = errorCovariance, verbose = TRUE, output_file_name.cis = output_file_name_cis, 
        pvOutputThreshold.cis = pvOutputThreshold_cis, snpspos = snpspos, 
        genepos = genepos, cisDist = cisDist, pvalue.hist = "qqplot", 
        min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)
    cat("Analysis done in: ", me$time.in.sec, " seconds", "\n")
    cat("Detected ", me$cis$neqtls, " local eQTLs:", "\n")
    cat("Detected ", me$trans$neqtls, " distant eQTLs:", "\n")
    plot(me)
    return(me)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
