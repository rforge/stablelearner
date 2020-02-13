%FIXME: is this the right way to document the S3 dispatch for RandomForest?

\name{stabletree-coercion}
\alias{stabletree-coercion}
\alias{as.stabletree}
\alias{as.stabletree.randomForest}
\alias{as.stabletree.RandomForest}
\alias{as.stabletree,RandomForest-method}
\alias{as.stabletree.cforest}

\title{Coercion Functions}

\description{
  Functions coercing various forest objects to objects of class
  \code{"stabletree"}.
}

\usage{
  as.stabletree(x, \dots)
  \method{as.stabletree}{randomForest}(x, applyfun = NULL, cores = NULL, \dots)
  \method{as.stabletree}{RandomForest}(x, applyfun = NULL, cores = NULL, \dots)
  \method{as.stabletree}{cforest}(x, applyfun = NULL, cores = NULL, savetrees = FALSE, \dots)
}

\arguments{
  \item{x}{an object of class \code{\link[randomForest]{randomForest}},
    \code{\link[party]{RandomForest-class}} or \code{\link[partykit]{cforest}}.}
  \item{applyfun}{a \code{\link[base]{lapply}}-like function. The default is to
    use \code{\link[base]{lapply}} unless \code{cores} is specified in which
    case \code{\link[parallel]{mclapply}} is used (for multicore computations on
    platforms that support these).}
  \item{cores}{integer. The number of cores to use in multicore computations
    using \code{\link[parallel]{mclapply}} (see above).}
  \item{savetrees}{logical. If \code{TRUE}, trees of the forests are returned.}
  \item{\dots}{additional arguments (currently not used).}
}

\details{
  Random forests fitted using \code{\link[randomForest]{randomForest}},
  \code{\link[party]{cforest}} or \code{\link[partykit]{cforest}} are coerced to
  \code{"stabletree"} objects.

  Note that when plotting a \code{\link[randomForest]{randomForest}}, the gray
  areas of levels of a nominal variable do not mimic exactly the same
  behavior as for classical \code{"stabletree"} objects, due to
  \code{\link[randomForest]{randomForest}} not storing any information whether
  any individuals were left fulfilling the splitting criterion in the subsample.
  Therefore, gray areas only indicate that this level of this variable has
  already been used in a split before in such a way that it could not be used
  for any further splits.
}

\value{
  \code{as.stabletree} returns an object of class \code{"stabletree"} which is a
  list with the following components:
  \item{call}{the call from the model object \code{x},}
  \item{B}{the number of trees of the random forest,}
  \item{sampler}{the random forest fitting function,}
  \item{vs0}{FIXME better name?,} 
  \item{br0}{FIXME better name?,}
  \item{vs}{FIXME better name?,}
  \item{br}{FIXME better name?,}
  \item{classes}{character vector indicating the classes of all partitioning
    variables,}
  \item{trees}{a list of tree objects of class \code{"party"}, or \code{NULL}.}
}

\seealso{\code{\link{stabletree}}, \code{\link{plot.stabletree}}}

\examples{

\donttest{
## build a randomForest using randomForest
library(randomForest)
set.seed(1)
rf <- randomForest(Species ~ ., data = iris)

## coerce to a stabletree
srf <- as.stabletree(rf)
print(srf)
summary(srf, original = FALSE) # there is no original tree
barplot(srf)
image(srf)
plot(srf)

## build a RandomForest using party
library(party)
set.seed(2)
cf_party <- cforest(Species ~ ., data = iris,
  control = cforest_unbiased(mtry = 2))

## coerce to a stabletree
scf_party <- as.stabletree(cf_party)
print(scf_party)
summary(scf_party, original = FALSE)
barplot(scf_party)
image(scf_party)
plot(scf_party)

## build a cforest using partykit
library(partykit)
set.seed(3)
cf_partykit <- cforest(Species ~ ., data = iris)

## coerce to a stabletree
scf_partykit <- as.stabletree(cf_partykit)
print(scf_partykit)
summary(scf_partykit, original = FALSE)
barplot(scf_partykit)
image(scf_partykit)
plot(scf_partykit)
}

}

\keyword{regression}

