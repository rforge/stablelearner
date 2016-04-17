\name{bootstrap}
\alias{bootstrap}
\alias{crossvalidation}
\alias{jackknife}
\alias{subsampling}

\title{Sampler Infrastructure for Stability Assessment}

\description{
  FIXME.
}

\usage{
  bootstrap(B = 500)
  subsampling(B = 500, frac = 0.8)
  crossvalidation(k = 5)
  jackknife(d = 1, maxrep = 5000)
}

\arguments{
  \item{B}{integer. Number of resampled datasets.}
  \item{k}{integer. Number of cross-validation folds.}
  \item{d}{integer. Number of FIXME.}
  \item{maxrep}{integer. Number of FIXME.}
  \item{frac}{numeric. Fraction of observations in
    each subsample.}
}

\details{
  FIXME
}

\seealso{\code{\link{stabletree}}}

\examples{
set.seed(0)

## bootstrap sampler
s <- bootstrap(3)
s$sampler(10)

## subsampling
s <- subsampling(3, frac = 0.6)
s$sampler(10)

## 5-fold cross-validation
s <- crossvalidation(5)
s$sampler(10)

## jackknife
s <- jackknife(d = 1)
s$sampler(10)
}

\keyword{regression}
