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
mxeqtl(snp_file, snp_location, expr_file, expr_location, 
cis_output_file, cis_pval, covariates = "", trans_output_file = "",
trans_pval = 0, model = "linear", MAF = 0, cis_dist = 1e+06,
missing="NA")
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
  \item{missing}{
    String for missing data. Default is NA.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
