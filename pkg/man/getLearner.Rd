\name{getLearner}
\alias{getLearner}

\encoding{latin1}

\title{Get learner details from \code{\link{LearnerList}}}

\description{
  Function to get information available about a specific learner in 
  \code{\link{LearnerList}} of the current \R session.
}

\usage{getLearner(x)}

\arguments{
  \item{x}{a fitted model object.}
}

\details{
  The function returns the entry in \code{\link{LearnerList}} found for 
  the class of the object submitted to the function.
}

\seealso{\code{\link{LearnerList}}, \code{\link{addLearner}}}

\examples{
m <- ctree(Species ~ ., data = iris)
getLearner(m)
}

\keyword{stability, similarity, resampling}
