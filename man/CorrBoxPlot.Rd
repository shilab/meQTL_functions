\name{CorrBoxPlot}
\alias{CorrBoxPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CorrBoxPlot
%%  ~~function to do ... ~~
}
\description{Creates Box Plots
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrBoxPlot(mEQTL, threshold, expr, genot, visual = FALSE, pdf_file = "", cis = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mEQTL}{
%%     ~~Describe \code{mEQTL} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{expr}{
%%     ~~Describe \code{expr} here~~
}
  \item{genot}{
%%     ~~Describe \code{genot} here~~
}
  \item{visual}{
%%     ~~Describe \code{visual} here~~
}
  \item{pdf_file}{
%%     ~~Describe \code{pdf_file} here~~
}
  \item{cis}{
%%     ~~Describe \code{cis} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
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
function (mEQTL, threshold, expr, genot, visual = FALSE, pdf_file = "", 
    cis = TRUE) 
{
    corr <- NULL
    phenotype <- NULL
    genotype <- NULL
    if (cis == TRUE) {
        index <- which(mEQTL$cis$eqtls$FDR <= threshold)
        eqtls <- mEQTL$cis$eqtls[index, ]
    }
    else {
        index <- which(mEQTL$trans$eqtls$FDR <= threshold)
        eqtls <- mEQTL$trans$eqtls[index, ]
    }
    for (i in 1:nrow(eqtls)) {
        phenotype[[i]] <- expr[which(expr$id == as.character(eqtls$gene[i])), 
            2:ncol(expr)]
        genotype[[i]] <- genot[which(genot$id == as.character(eqtls$snps[i])), 
            2:ncol(genot)]
        corr[i] <- cor(as.numeric(phenotype[[i]]), as.numeric(genotype[[i]]), 
            use = "pairwise.complete.obs")
    }
    if (visual) {
        pdf(paste("./", pdf_file, sep = ""))
        par(mfcol = c(2, 2))
        genotypes <- range(genot[, 2:ncol(genot)])[1]:range(genot[, 
            2:ncol(genot)])[2]
        for (i in 1:nrow(eqtls)) {
            pheno <- as.numeric(phenotype[[i]])
            values <- NULL
            for (j in 1:length(genotypes)) {
                values[[j]] <- pheno[which(genotype[[i]] == genotypes[j])]
            }
            if (length(genotypes) == 3) {
                cats = c(0, 1, 2)
            }
            else {
                cats = c(0, 1, 2, 3, 4)
            }
            boxplot(values, boxwex = 0.5, ylab = paste(as.character(eqtls$gene[i]), 
                " expression"), names = cats, xlab = paste(as.character(eqtls$snps[i]), 
                " genotype", "\nCorrelation: ", format(corr[i], 
                  2), "P-value: ", format(eqtls$pvalue[i], 2), 
                " FDR: ", format(eqtls$FDR[i], 2)), main = paste(as.character(eqtls$snps[i]), 
                " - ", as.character(eqtls$gene[i])))
        }
        dev.off()
    }
    return(corr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
