###Installing
meQTLfunc can be installed by using the devtools package.
```
install_github('shilab/meQTL_functions')
```

``` {.r}
library(MatrixEQTL)
library(meQTLfunc)
```

``` {.r}
setwd('~/Development/Repos/meQTLfunctions/')
me<-mxeqtl('CNV_matrix.out.filter','CNV_position','liver_expression.out.filter','gene_position','liver_cis_results',0.05)
```

    ## Rows read:  2,000 
    ## Rows read:  2529  done.
    ## SNPs before filtering: 2529 
    ## Rows read:  2,000 
    ## Rows read:  4,000 
    ## Rows read:  4557  done.
    ## Matching data files and location files 
    ## 2591 of 4557  genes matched
    ## 2529 of 2529  SNPs matched
    ## Task finished in  0.353  seconds
    ## Reordering SNPs
    ##  
    ## Task finished in  0.094  seconds
    ## Reordering genes
    ##  
    ## Task finished in  0.098  seconds
    ## Processing covariates 
    ## Task finished in  0.003  seconds
    ## Processing gene expression data (imputation, residualization, etc.) 
    ## Task finished in  0.017  seconds
    ## Creating output file(s) 
    ## Task finished in  0.009  seconds
    ## Performing eQTL analysis 
    ## 16.66% done, 585 cis-eQTLs
    ## 66.66% done, 605 cis-eQTLs
    ## 83.33% done, 780 cis-eQTLs
    ## Task finished in  0.91  seconds
    ##  
    ## Analysis done in:  1.457  seconds 
    ## Detected  780  local eQTLs: 
    ## Detected   distant eQTLs:

![plot of chunk unnamed-chunk-2](./README_files/figure-markdown_github/unnamed-chunk-2.png)

``` {.r}
expr<-read.table('./liver_expression.out.filter', header = TRUE, stringsAsFactors = FALSE,na.string="NA")
genot<-read.table('./CNV_matrix.out.filter', header = TRUE, stringsAsFactors = FALSE,na.string="NA")
CorrBoxPlot(me,.2,expr,genot,visual=T,pdf_file="res.pdf")
```

    ##  [1]  0.7859  0.7799  0.7799 -0.7799  0.7554 -0.7412 -0.7301  0.7113
    ##  [9]  0.7015  0.6977 -0.6936  0.6935  0.6902 -0.6845  0.6840 -0.6839
    ## [17]  0.6828  0.6782  0.6773  0.6745  0.6744  0.6744 -0.6710  0.6701
    ## [25]  0.6693  0.6685  0.6658  0.6634  0.6608 -0.6540  0.6524 -0.6524
    ## [33]  0.6524 -0.6524  0.6480 -0.6463  0.6422  0.6403  0.6393  0.6383
    ## [41] -0.6345  0.6345 -0.6345 -0.6345 -0.6336  0.6333  0.6333 -0.6333
    ## [49]  0.6328  0.6306  0.6289 -0.6258  0.6256 -0.6256 -0.6256  0.6246
    ## [57] -0.6228  0.6226  0.6221  0.6196  0.6188  0.6183  0.6183 -0.6183
    ## [65]  0.6155  0.6148  0.6143 -0.6118 -0.6099 -0.6095  0.6093  0.6083
    ## [73] -0.6078  0.6075  0.6065  0.6063

``` {.r}
CorrBoxPlot
```

    ## function (mEQTL,threshold,expr,genot,visual=FALSE,pdf_file="",crlt=0,cis=TRUE){
    ##     # Inputs:
    ##     #   mEQTL     - Matrix EQTL object with the eQTLs already collected
    ##     #   threshold - FDR cutoff, only those eQTLs with equal or lower threshold will be taken into account
    ##     #   expr      - Transcript expression dataset
    ##     #   genot     - Genotyping dataset, either phased or unphased variants
    ##     #   cis       - If TRUE only cis eQTLs are considered, otherwise trans eQTLS (default TRUE)
    ##     #   visual    - If TRUE the script will display a box plot figure for each eQTL above the threshold
    ##     #
    ##     # Output: A vector with Pearson correlation scores for each eQTL that surpasses the given threshold
    ##     #
    ##     # NOTES:
    ##     # Obviously the original files from which mEQTL object was computed must match on transcript, variants
    ##     # and samples IDs included in expr and genot 
    ##     #
    ##     # expr and genot datafiles are in the matrixEQTL format and can be loaded as:
    ##     # expr = read.table(file_name, header = TRUE, stringsAsFactors = FALSE);
    ##     #
    ##     # R. Armananzas - Last update 12/02/13
    ##     #
    ##     
    ##     
    ##     corr  <- NULL; phenotype <- NULL; genotype <- NULL; phenotype_new<-NULL;
    ##     
    ##     if (cis==TRUE){
    ##       index <- which(mEQTL$cis$eqtls$FDR<=threshold)
    ##       eqtls <- mEQTL$cis$eqtls[index,]
    ##     }
    ##     else{
    ##       index <- which(mEQTL$trans$eqtls$FDR<=threshold)
    ##       eqtls <- mEQTL$trans$eqtls[index,]    
    ##     }
    ##     
    ##     for (i in 1:nrow(eqtls)){
    ##       phenotype[[i]] <- expr[which(expr[,1]==as.character(eqtls$gene[i])),2:ncol(expr)]
    ##       genotype[[i]]  <- genot[which(genot[,1]==as.character(eqtls$snps[i])),2:ncol(genot)]
    ##       corr[i]   <- cor(as.numeric(phenotype[[i]]),as.numeric(genotype[[i]]), use="pairwise.complete.obs")
    ##     }
    ##     
    ##     if (visual){ #Perform the plots
    ##       #There can be three values (unphased) or four (phased)
    ##       pdf(paste('./',pdf_file,sep=""))
    ##       par(mfcol = c(2, 2))
    ##       genotypes <- range(genot[,2:ncol(genot)])[1]:range(genot[,2:ncol(genot)])[2]
    ##       for (i in 1:nrow(eqtls)){
    ##         #Prepare the matrix
    ##         pheno <- as.numeric(phenotype[[i]])
    ##         values <- NULL;
    ##         for (j in 1:length(genotypes)){
    ##           values[[j]] <- pheno[which(genotype[[i]]==genotypes[j])]
    ##         }
    ##         #Plot the boxplots
    ##         #if (length(genotypes)==3){cats=c(0,1,2)}
    ##         #else {cats=c(0,1,2,3,4)}
    ##         cats=seq(0,length(genotypes)-1)
    ##         if (abs(corr[i])>=crlt)
    ##         {
    ##           boxplot(values,boxwex=0.5,ylab=paste(as.character(eqtls$gene[i])," expression"), names=cats,
    ##                   xlab=paste(as.character(eqtls$snps[i])," genotype","\nCorrelation: ",format(corr[i],2),
    ##                              "P-value: ",format(eqtls$pvalue[i],2)," FDR: ",format(eqtls$FDR[i],2)),
    ##                   main=paste(as.character(eqtls$snps[i])," - ",as.character(eqtls$gene[i])))
    ##         }
    ##       }
    ##       dev.off()
    ##     }
    ##     return(corr)
    ##   }
    ## <environment: namespace:meQTLfunc>
