\name{Bipolar2009}
\alias{Bipolar2009}

\title{Case-Control Study on Bipolar Disorder (2009 Data)}

\description{
  Data originally presented by Ryan, Lockstone, Huffaker, Wayland, Webster &
  Bahn (2006) from the dorsolateral prefrontal cortex of 61 persons (31
  controls) including their age in years, their brain pH-level, their status
  (bipolar disorder vs. healthy control), their gender (female vs. male). The
  modified dataset presented here was introduced by Strobl, Malley & Tutz (2009)
  and includes 100 randomly selected genes of the 22283 original genes as well
  as two genes that were artificially modified to have an effect on the status
  of the persons.
}

\usage{data("Bipolar2009", package = "stablelearner")}

\format{
A data frame containing 61 observations on 106 variables.
\describe{
  The data frame contains the variables \code{age} (age in years),
  \code{brain_pH} (numeric value of the brain pH-level), \code{status} (a
  factor with levels \code{Bipolar disorder} and \code{Healthy control}),
  \code{gender} (a factor with levels \code{Female} and \code{Male}) and 100
  randomly selected genes of the original dataset (prefixed with \code{gene_}),
  as well as two artificially modified genes to have an effect on the
  \code{status}, namely \code{simulated_gene_1} and \code{simulated_gene_2}
  (both of numeric value).
}
}

\details{
Ryan et al. (2006) used microarray technology to determine the expression of
approximately 22000 mRNA transcripts in post-mortem dorsolateral prefrontal
cortex tissue in patients with bipolar disorder (30 persons) and matched
healthy controls (31 persons). In the original study of Ryan et al. (2006),
no genes were clearly found to be differential expressed. To demonstrate the
application of random forests on this dataset, Strobl et al. (2009) randomly
selected 100 genes and modified two genes artificially to have an effect on the
status of the persons.
}

\references{
  Ryan MM, Lockstone HE, Huffaker SJ, Wayland MT, Webster MJ, Bahn S (2006).
    Gene Expression Analysis of Bipolar Disorder Reveals Downregulation of the
    Ubiquitin Cycle and Alterations in Synaptic Genes.
    \emph{Molecular psychiatry}, \bold{11}(10), 965--978.

  Strobl C, Malley J, Tutz G (2009).
    An Introduction to Recursive Partitioning: Rationale, Application, and
    Characteristics of Classification and Regression Trees, Bagging, and Random
    Forests.
    \emph{Psychological Methods}, \bold{14}(4), 323--348.
}

\examples{
## overview
data("Bipolar2009", package = "stablelearner")
str(Bipolar2009)
}

\keyword{datasets}

