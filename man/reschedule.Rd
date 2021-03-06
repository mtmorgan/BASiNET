% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reschedule.R
\name{reschedule}
\alias{reschedule}
\title{Rescales the results between values from 0 to 1}
\usage{
reschedule(matrix, mRNA, lncRNA, sncRNA)
}
\arguments{
\item{matrix}{Array with results numerics}

\item{mRNA}{Integer number of mRNA sequences}

\item{lncRNA}{Integer number of lncRNA sequences}

\item{sncRNA}{Integer number of sncRNA sequences}
}
\value{
Returns the array with the rescaled values
}
\description{
Given the results the data is rescaled for values between 0 and 1, so that the length of the sequences does not influence the results. The rescaling of the mRNA and lncRNA are made separately
}
\author{
Eric Augusto Ito
}
