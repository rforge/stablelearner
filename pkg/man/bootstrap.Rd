\name{bootstrap}
\alias{bootstrap}
\alias{crossvalidation}
\alias{jackknife}
\alias{subsampling}

\title{Sampler Infrastructure for Stability Assessment}

\description{
  Sampler objects that provide objects with functionality used by 
  \code{\link{stabletree}} to generate resampled datasets.
}

\usage{
  bootstrap(B = 500)
  subsampling(B = 500, frac = 0.8)
  crossvalidation(k = 5)
  jackknife(d = 1, maxrep = 5000)
}

\arguments{
  \item{B}{An integer value specifying the number of resampled datasets.}
  \item{k}{An integer value specifying the number of folds in cross validation.}
  \item{d}{An integer value specifying the number of observations left out in jackknife.}
  \item{maxrep}{An integer value specifying the maximum number of resampled datasets allowed, when using jackknife.}
  \item{frac}{A numeric value between 0 and 1 specifying the fraction of observations in each subsample.}
}

\details{

  The sampler functions provide objects that include functionality to generate 
  resampled datasets used by \code{\link{stabletree}}.

  The \code{bootstrap} function provides an object that can be used to generate
  \code{B} bootstrap samples by sampling from \code{n} observations with 
  replacement.

  The \code{subsampling} function provides an object that can be used to 
  generate \code{B} subsamples by sampling from \code{floor(frac*n)} 
  observations without replacement.

  The \code{crossvalidation} function provides an object that can be used to 
  generate \code{k}-folds from \code{n} observations.

  The \code{jackknife} function provides an object that can be used to generate
  all datasets necessary to perform leave-\code{k}-out jackknife sampling from 
  \code{n} observations. The number of datasets is limited by \code{maxrep} to 
  prevent unintended CPU or memory overload by accidently choosing too large 
  values for \code{k}.
  
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
