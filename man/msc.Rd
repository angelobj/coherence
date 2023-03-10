% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coherence.R
\name{msc}
\alias{msc}
\title{Phase angle between complex valued vectors (based on Dr. Christopher Laine matlab code).}
\usage{
msc(x, y)
}
\arguments{
\item{x}{complex valued vector}

\item{y}{complex valued vector}
}
\description{
Phase angle between complex valued vectors (based on Dr. Christopher Laine matlab code).
}
\examples{
x<-c(2 +   2i,3 +   2i,4 +   7i,5 +   7i)
y<-c(3 +   2i,4 +   2i,5 +   7i,6 +   7i)
msc(x,y)
}
