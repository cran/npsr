% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npsValid.R
\name{nps.valid}
\alias{nps.valid}
\title{nps.valid}
\usage{
nps.valid(Q, l, m, n, N = sum(Q), S = sum(Q))
}
\arguments{
\item{Q}{Histogram of dataset (l*m*n vector)}

\item{l}{|Z|}

\item{m}{|X|}

\item{n}{|Y|}

\item{N}{Number of Repetitions for Nested Sampling}

\item{S}{Number of Starting Points for Nested Sampling}
}
\description{
Calculates M_Valid
}
