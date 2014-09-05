\name{CorrScatterPlot}
\alias{CorrScatterPlot}
%- Also NEED an '\alias' for EACH other topic documented here
\title{ Correlation Scatter Plots
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrScatterPlot(mEQTL, threshold, expr, genot, visual = TRUE, cis = TRUE)
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
function (mEQTL, threshold, expr, genot, visual = TRUE, cis = TRUE) 
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
        phenotype[[i]] <- expr[which(expr$geneid == as.character(eqtls$gene[i])), 
            2:ncol(expr)]
        genotype[[i]] <- genot[which(genot$snpid == as.character(eqtls$snps[i])), 
            2:ncol(genot)]
        corr[i] <- cor(as.numeric(phenotype[[i]]), as.numeric(genotype[[i]]))
    }
    if (visual) {
        genotypes <- range(genot[, 2:ncol(genot)])[1]:range(genot[, 
            2:ncol(genot)])[2]
        for (i in 1:nrow(eqtls)) {
            geno <- as.numeric(genotype[[i]])
            pheno <- as.numeric(phenotype[[i]])
            values <- NULL
            for (j in 1:length(genotypes)) {
                values[[j]] <- pheno[which(genotype[[i]] == genotypes[j])]
            }
            plot(values, xlab = "miRNA Expression", ylab = "mRNA Expression", 
                main = paste(as.character(eqtls$snps[i]), " genotype", 
                  "\nCorrelation: ", format(corr[i], 2), "P-value: ", 
                  format(eqtls$pvalue[i], 2), " FDR: ", format(eqtls$FDR[i], 
                    2)))
        }
    }
    return(corr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
