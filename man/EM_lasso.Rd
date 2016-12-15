\name{EM_lasso}
\alias{EM_lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The EM algorithm for the Spike and Slab model inference
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EM_lasso(S, n, p_n, v0, v1, maxiter, p, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
Sample Covariance Matrix
}
  \item{n}{
%%     ~~Describe \code{n} here~~
Number of observations
}
  \item{p_n}{
%%     ~~Describe \code{p_n} here~~
Number of variabes/ncol in the precision matrix
}
  \item{v0}{
%%     ~~Describe \code{v0} here~~
The spike tuning parameter
}
  \item{v1}{
%%     ~~Describe \code{v1} here~~
The slab tuning parameter
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
Max iteration ties
}
  \item{p}{
%%     ~~Describe \code{p} here~~
eta in the paper: the prior probability of non sparsity
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
The penalty on the diagnoal term
}
}

\value{
	\item{P}{Marginal posterior probabilty matrix of each entry being sparse}
	\item{Theta}{An MAP estimate of the precision matrix}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
