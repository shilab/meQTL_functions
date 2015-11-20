[![Build Status](https://travis-ci.org/shilab/meQTL_functions.svg?branch=master)](https://travis-ci.org/shilab/meQTL_functions) [![Coverage Status](https://coveralls.io/repos/shilab/meQTL_functions/badge.svg?branch=master&service=github)](https://coveralls.io/github/shilab/meQTL_functions?branch=master)

###Installing
meQTLfunc can be installed by using the devtools package.
```
install_github('shilab/meQTL_functions')
library(meQTLfunc)
```

``` {.r}
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


