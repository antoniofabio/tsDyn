\name{autotriples.rgl}
\alias{autotriples.rgl}
\title{Interactive trivariate time series plots}
\description{
Interactive trivariate time series plots
}
\usage{
autotriples.rgl(x, lags=1:2, type=c("lines", "points"))
}
\arguments{
\item{x}{ time series}
\item{lags}{ vector of regressors lags }
\item{type}{ type of plot: countour levels, perspective plots, image }
}
\details{
This function displays interactive trivariate time series plots 
\code{x[t-lags[1]], x[t-lags[2]]} against \code{x[t]} using the interactive \code{\link[rgl]{rgl}} device.
}
\value{
None. A plot is produced on the current \code{rgl} device.
}
\seealso{
\code{\link{autotriples}} for 3d visualization via \code{scatterplot3d} package and for kernel post-processing of the cloud for nonparametric autoregression functions estimates.
}
\examples{
if(interactive())
 autotriples.rgl(log(lynx))
}
\author{ Wrapper to 'sm' and GUI by Antonio, Fabio Di Narzo  }
\keyword{ ts }
