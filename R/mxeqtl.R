toSkip <- function(exists)
{
    if (exists)
    {
        return(1);
    }
    else
    {
        return(0);
    }
}

getGenotypes <- function(sep, missing, header, rownames, snp_filename)
{
    if (!file.exists(snp_filename))
    {
        stop(cat(snp_filename,' does not exist','\n'));
    }
    skipRows <- toSkip(header);
    skipCols <- toSkip(rownames);
    snps <- SlicedData$new();
    snps$fileDelimiter <- sep;
    snps$fileOmitCharacters <- missing;
    snps$fileSkipRows <- skipRows;
    snps$fileSkipColumns <- skipCols;
    snps$fileSliceSize <- 2000;
    snps$LoadFile(snp_filename);
    return(snps);
}

getExpression <- function(sep, missing, header, rownames, expr_filename)
{
    if (!file.exists(expr_filename))
    {
        stop(cat(expr_filename,' does not exist','\n'));
    }
    skipRows <- toSkip(header);
    skipCols <- toSkip(rownames);    
    gene <- SlicedData$new();
    gene$fileDelimiter <- sep;
    gene$fileOmitCharacters <- missing;
    gene$fileSkipRows <- skipRows;
    gene$fileSkipColumns <- skipCols;
    gene$fileSliceSize <- 2000;
    gene$LoadFile(expr_filename);
    return(gene)
}

getCovariates <- function(sep, missing, header, rownames, cvrt_filename)
{
    skipRows <- toSkip(header);
    skipCols <- toSkip(rownames);
    cvrt <- SlicedData$new();
    cvrt$fileDelimiter <- sep;
    cvrt$fileOmitCharacters <- missing;
    cvrt$fileSkipRows <- skipRows;
    if (length(cvrt_filename) == 0)
    {
        return(cvrt);
    }
    if (!file.exists(cvrt_filename))
    {
        stop(cat(cvrt_filename,' does not exist','\n'));
    }
    cvrt$LoadFile(cvrt_filename);
    return(cvrt)
}

mafFilter <- function(snps, MAF)
{
    Mode <- function(x) 
    {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }

    cat('SNPs before filtering:',nrow(snps),'\n')
    maf.list <- vector('list', length(snps))
    for(sl in 1:length(snps)) 
    {
        slice <- snps[[sl]];
        maf.list[[sl]] <- rowMeans(slice != apply(slice, 1, Mode),na.rm=T)/2;
    }
    maf <- unlist(maf.list)
    sum(maf>=MAF)
    snps$RowReorder(maf>MAF);
    cat('SNPs after filtering:',nrow(snps),'\n');
    rm(maf, sl, maf.list);    
    return(snps)
}

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

getCovFileName <-function(covariates)
{
    if (covariates!="")
    {
        return(covariates);
    }
    else
    {
        return(character());
    }
}

#' @export
mxeqtl <-
function(snp_file,snp_location,expr_file,expr_location,cis_output_file,
         cis_pval,covariates="",trans_output_file="", trans_pval=0, 
         model="linear", MAF=0, cis_dist=1e6, qq="", missing="NA", sep="\t", 
         header=TRUE,rownames=TRUE)
{
    # Matrix eQTL function based on the sample code by Andrey A. Shabalin
    # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

    useModel <- setModel(model) #nocov start

    # Covariates file name

    covariates_file_name <- getCovFileName(covariates);

    # Output file name
    output_file_name_cis <- cis_output_file
    if (trans_output_file!="")
    {
        output_file_name_tra <- trans_output_file;
    }
    else
    {
        output_file_name_tra <- "onlyTRANSresults.txt"
    }

    #TODO: Create option for errorCovariance
    errorCovariance <- numeric();

    snps <- getGenotypes(sep, missing, header, rownames, snp_file);

    if (MAF>0)
    {
        snps <- mafFilter(snps, MAF);
	    snps$SaveFile("meQTL_filtered_input")
    }

    gene <- getExpression(sep, missing, header, rownames, expr_file);

    cvrt <- getCovariates(sep, missing, header, rownames, covariates_file_name);

    # Load the genotype and expression positions
    snpspos <- read.table(snp_location, header = TRUE, stringsAsFactors = FALSE);
    genepos <- read.table(expr_location, header = TRUE, stringsAsFactors = FALSE);

    #Call MatrixEQTL function
    me <- Matrix_eQTL_main(
    snps <- snps,
    gene <- gene,
    cvrt <- cvrt,
    output_file_name     <- output_file_name_tra,
    pvOutputThreshold     <- trans_pval,
    useModel <- useModel,
    errorCovariance <- errorCovariance,
    verbose <- TRUE,
    output_file_name.cis <- output_file_name_cis,
    pvOutputThreshold.cis <- cis_pval,
    snpspos <- snpspos,
    genepos <- genepos,
    cisDist <- cis_dist,
    pvalue.hist <- "qqplot",
    min.pv.by.genesnp <- TRUE,
    noFDRsaveMemory  <- FALSE);

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

    return(me) #nocov end
}
