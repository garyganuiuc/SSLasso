\name{Tune_SSLasso}
\alias{Tune_SSLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Tuning SSLasso
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Tune_SSLasso(v0, tau, S, n, p_n, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v0}{
  	Tuning parameter v0
%%     ~~Describe \code{v0} here~~
}
  \item{tau}{
  	  	Tuning parameter tau, control on the diagonal
%%     ~~Describe \code{tau} here~~
}
  \item{S}{
  	Sample covariance matrix
%%     ~~Describe \code{S} here~~
}
  \item{n}{
  	Sample size
%%     ~~Describe \code{n} here~~
}
  \item{p_n}{
  	Matrix dimension
%%     ~~Describe \code{p_n} here~~
}
  \item{p}{
  	The hyperparameter eta
%%     ~~Describe \code{p} here~~
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
function (v0, tau, S, n, p_n, p) 
{
    pb <- txtProgressBar(min = 0, max = length(v0) * 10 * length(tau), 
        style = 3)
    bic = NULL
    for (m in 1:length(p)) {
        for (i in 1:length(v0)) {
            v1 = seq(v0[i] + 0.1, 5 * v0[i], 0.5 * v0[i])
            for (j in 1:length(v1)) {
                for (k in 1:length(tau)) {
                  w = 1
                  l = 1
                  maxiter = 30
                  result1 <- EM_lasso(S, n, p_n, v0[i], v1[j], 
                    maxiter, p[m], tau[k])
                  bic = rbind(bic, list(v0 = v0[i], v1 = v1[j], 
                    tau = tau[k], p = p[m], BIC = BIC_SSLasso(result1$Theta, 
                      S, result1$P, n)))
                  setTxtProgressBar(pb, (m - 1) * length(v0) * 
                    length(v1) * length(tau) + (i - 1) * length(v1) * 
                    length(tau) + (j - 1) * length(tau) + k)
                }
            }
        }
    }
    close(pb)
    return(bic[which.min(bic[, 5]), ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
