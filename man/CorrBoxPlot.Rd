\name{CorrBoxPlot}
\alias{CorrBoxPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates Correlation Box Plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrBoxPlot(mEQTL, threshold, expr, genot, visual = FALSE,
pdf_file = "", crlt = 0, cis = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mEQTL}{
The name of the variable containing the results from mxeqtl. Required.
}
  \item{threshold}{
FDR cutoff. Required.
}
  \item{expr}{
The file containing the expression data. Needs to be loaded:
expr = read.table(file_name, header = TRUE, stringsAsFactors = FALSE,na.string="NA");
Required.
}
  \item{genot}{
The file containing the genotype data. Needs to be loaded:
genot = read.table(file_name, header = TRUE, stringsAsFactors = FALSE,na.string="NA");
Required.
}
  \item{visual}{
Set to TRUE to print plots to PDF file. Optional.
}
  \item{pdf_file}{
Name of PDF file for plots. Optional
}
  \item{crlt}{
Correlation cutoff, for filtering plots. Default is 0, no filtering. Optional
}
  \item{cis}{
Set to FALSE if results contain trans eQTLs. 
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
