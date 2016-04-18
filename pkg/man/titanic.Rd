\name{titanic}
\alias{titanic}

\title{Passengers and Crew on the RMS Titanic}

\description{
  the \code{titanic} data is a complete list of passengers and crew members on 
  the RMS Titanic. It includes a variable indicating whether a person did 
  survive the sinking of the RMS Titanic on 15th April 1912.
}

\usage{data("titanic")}

\format{
A data frame containing 2207 observations on 11 variables.
\describe{
  \item{name}{.} a character with the persons full name.
  \item{gender}{.} a factor with levels \code{male} and \code{female}.
  \item{age}{.} a numeric value with the age of the person on the day of the 
  sinking. The age of babies (under 12 months) is given as a fraction of one 
  year (1/month).
  \item{class}{.} a factor with the class for passengers (\code{1st}, 
  \code{2nd}, \code{3rd}) and the type of service aboard for crew members 
  (\code{deck crew}, \code{engineering crew}, \code{restaurant staff}, 
  \code{victualling crew}).
  \item{embarked}{.} a factor with the place of of embarkment 
  (\code{Belfast}, \code{Cherbourg}, \code{Queenstown}, \code{Southampton}).
  \item{country}{.} a factor with the persons home country.
  \item{ticketno}{.} a numeric value specifying the ticket number 
  (\code{NA} for crew members).
  \item{fare}{.} a numeric value with the ticket price (\code{NA} for crew 
  members, musicians and employees of the shipyard company).
  \item{sibsp}{.} an ordered factor specifying the number if siblings/spouses 
  aboard (adopted from Vanderbild data set).
  \item{parch}{.} an ordered factor specifying the number if parents/children 
  aboard (adopted from Vanderbild data set).
  \item{survived}{.} a factor with two levels (\code{no}, \code{yes}) 
  specifying whether the person has survived the sinking.
}
}

\details{

The website \url{http://www.encyclopedia-titanica.org/} offers detailed 
information about passengers and crew members on the RMS Titanic. According to
the website $1317$ passengers and $890$ crew member were abord.

$8$ musicians and $9$ employees of the shipyard company are listed as 
passengers, but travelled with a free ticket, which is why they have \code{NA}
values in \code{fare}. A few regular passengers also have missing values.

}

\source{
  The complete list of persons on the RMS titanic was downloaded from 
  \url{http://www.encyclopedia-titanica.org/} on April 5th, 2016. The 
  information given in \code{sibsp} and \code{parch} was adopoted from a data
  set obtained from \urk{http://biostat.mc.vanderbilt.edu/DataSets}.
}

\references{
  \url{http://www.encyclopedia-titanica.org/} and 
  \urk{http://biostat.mc.vanderbilt.edu/DataSets}.
}

\examples{
data("titanic", package = "stablelearner")
summary(titanic)

## passengers only
titanic_passenger <- subset(titanic, class %in% c("1st", "2nd", "3rd"))
summary(titanic_passenger)

## crew members only
titanic_crew <- subset(titanic, !c(class %in% c("1st", "2nd", "3rd")))
summary(titanic_crew)

}

\keyword{datasets}
