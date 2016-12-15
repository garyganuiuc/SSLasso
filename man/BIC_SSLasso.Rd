\name{BIC_SSLasso}
\alias{BIC_SSLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIC for SSLasso model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the BIC for the model selected.
}
\usage{
BIC_SSLasso(Theta, S, P, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Theta}{
%%     ~~Describe \code{Theta} here~~
The precision matrix estimated
}
  \item{S}{
%%     ~~Describe \code{S} here~~
Sample covariance matrix
}
  \item{P}{
%%     ~~Describe \code{P} here~~
The sparsity pattern matrix
}
  \item{n}{
%%     ~~Describe \code{n} here~~
Number of observations
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
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
