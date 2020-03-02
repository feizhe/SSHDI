\name{SSHDI}
\alias{SSHDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SSHDI
}
\description{
Fit SSHDI for high dimensional estimation and inference}
\usage{
SSHDI(xmat,yvec,family = "gaussian",B=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmat}{ Design matrix (no intercept).
%%     ~~Describe \code{x} here~~
}
  \item{yvec}{ Response variable.
%%     ~~Describe \code{x} here~~
}
  \item{family}{ family from c("gaussian","poisson", "multinomial",
  "cox", "mgaussian").
%%     ~~Describe \code{x} here~~
}
  \item{B}{ Number of re-samples. Recommend \code{B} equals sample size, or similar.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{int }{Intercept estimator.}
  \item{ss.beta }{Coefficient vector.}
  \item{sd}{Standard errors of the coefficients.}
  \item{p}{Individual p-values for each predictor.}
  \item{sel.freq}{Selection frequency for each predictor based on \code{B} re-samples.}
  \item{nzero}{Number of significant predictors after Bonferroni correction.}
%% ...
}
\references{
https://arxiv.org/abs/1903.04408
%% ~put references to the literature/web site here ~
}
\author{
Zhe Fei
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
