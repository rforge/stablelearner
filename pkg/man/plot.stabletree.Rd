\name{plot.stabletree}
\alias{plot.stabletree}
\alias{barplot.stabletree}
\alias{image.stabletree}

\title{Visualizing Tree Stability Assessments}

\description{
  Visualizations of tree stability assessments carried out
  via \code{\link{stabletree}}.
}

\usage{
\method{plot}{stabletree}(x, select = order(colMeans(x$vs), decreasing = TRUE), 
  type.breaks = "levels", col.breaks = "red", lty.breaks = "dashed",
  cex.breaks = 0.7,  col.main = c("black", "gray50"), main.uline = TRUE,
  args.numeric = NULL, args.factor = NULL, args.ordered = NULL, \dots)

\method{barplot}{stabletree}(height, main = "Variable selection frequencies",
  xlab = "", ylab = "", horiz = FALSE, col = gray.colors(2),
  names.arg = NULL, names.uline = TRUE, names.diag = TRUE,
  cex.names = 0.9, line = 0.5,
  ylim = if (horiz) NULL else c(0, 100), xlim = if (horiz) c(0, 100) else NULL,
  \dots)

\method{image}{stabletree}(x, main = "Variable selections",
  ylab = "Repetitions", xlab = "", col = gray.colors(2),
  names.arg = NULL, names.uline = TRUE, names.diag = TRUE, 
  cex.names = 0.9, xaxs = "i", yaxs = "i",
  xlim = c(0, length(x$vs0)), ylim = c(0, x$B), \dots)
}

\arguments{
  \item{x, height}{an object of class \code{\link{stabletree}}.}
  \item{select}{integer or character.}
  \item{type.breaks}{.}
  \item{col.breaks}{.}
  \item{lty.breaks}{.}
  \item{cex.breaks}{.}
  \item{col.main}{.}
  \item{main.uline}{.}
  \item{args.numeric}{.}
  \item{args.factor}{.}
  \item{args.ordered}{.}
  \item{\dots}{further arguments passed to plotting functions,
    especially for labeling and annotation.}
  \item{main, xlab, ylab}{character. Annotations of axes and main
    title, respectively.}
  \item{horiz}{.}
  \item{col}{.}
  \item{names.arg}{.}
  \item{names.uline}{.}
  \item{names.diag}{.}
  \item{cex.names}{.}
  \item{line}{.}
  \item{xlim, ylim}{.}
  \item{xaxs, yaxs}{.}
}

\details{
  FIXME
}

\seealso{\code{\link{stabletree}}}

\examples{
## build a tree
m <- ctree(Species ~ ., data = iris)
plot(m)

## investigate stability
set.seed(0)
s <- stabletree(m, B = 500)

## show variable selection proportions
## with different labels and different ordering
barplot(s)
barplot(s, cex.names = 0.8)
barplot(s, names.diag = FALSE)
barplot(s, names.arg = c("a", "b", "c", "d"))
barplot(s, names.uline = FALSE)
barplot(s, col = c("lightgreen", "darkred"))
barplot(s, horiz = TRUE)

## illustrate variable selections of replications
## with different labels and different ordering
image(s)
image(s, cex.names = 0.8)
image(s, names.diag = FALSE)
image(s, names.arg = c("a", "b", "c", "d"))
image(s, names.uline = FALSE)
image(s, col = c("lightgreen", "darkred"))
image(s, order.row = FALSE)
image(s, order.col = FALSE)

## graphical cutpoint analysis, selecting variable by number and name
## with different numerical of break points
plot(s)
plot(s, select = 3)
plot(s, select = "Petal.Width")
plot(s, args.numeric = list(breaks = 40))

# change labels of splits in complete data tree
plot(s, select = 3, type.labels = "levels")
plot(s, select = 3, type.labels = "nodeids")
plot(s, select = 3, type.labels = "breaks")
plot(s, select = 3, type.labels = "none")
}

\keyword{regression}
