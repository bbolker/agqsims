% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/AGQfuns.R
\name{simfun_culc}
\alias{simfun_culc}
\title{simulation code, based approximately on Culcita data set
not longitudinal; approximates toenail if n.ttt=2
simulates covariate x, but not currently used}
\usage{
simfun_culc(n.blocks = 10, n.ttt = 4, n.rep = 3, N.binom = 1,
  x.range = c(-1, 1), beta = c(5, -3.75, -4.4, -5.5), theta = 3.5,
  seed = NULL)
}
\arguments{
\item{n.blocks}{number of levels of RE grouping variable}

\item{n.ttt}{number of levels of (categorical) fixed effect}

\item{n.rep}{number of replicates per ttt*block combination}

\item{N.binom}{number of trials per binomial sample}

\item{x.range}{range of continuous covariate (unused)}

\item{beta}{fixed-effects parameter vector}

\item{theta}{RE parameter vector (Cholesky factor)}

\item{seed}{random-number seed}
}
\description{
simulation code, based approximately on Culcita data set
not longitudinal; approximates toenail if n.ttt=2
simulates covariate x, but not currently used
}

