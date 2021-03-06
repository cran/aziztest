% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_analysis.R
\name{which_aberrant}
\alias{which_aberrant}
\title{Returns which samples are in the aberrant interval defined by \code{test.aziz()}}
\usage{
which_aberrant(xi, x, res)
}
\arguments{
\item{xi}{Numerical vector. Can be the same as the tested variable x
or it can be a new set of unseen samples.}

\item{x}{Numerical vector of the variable tested by \code{test.aziz()}}

\item{res}{Result of running \code{test.aziz()}}
}
\value{
indexes of samples in \code{xi} that are within the aberrant interval
}
\description{
Returns index of samples that are in the aberrant interval defined by \code{test.aziz()}
Can take the same samples or a new set of previously unseen samples
}
\examples{
y = c(rep(1,200),rep(0,200))
x = rnorm(400)
#Inducing an aberration enrichment signal by perturbing some of the cases
x[1:20]=x[1:20]-3;
res2 = aziz.test(y,x,rep=20000)
print_summary(res2)
which_aberrant(x,x,res2)
which_aberrant(c(-5,1.5,-2.5,-0.5,2),x,res2)#testing if new values are within the aberrant interval
}
