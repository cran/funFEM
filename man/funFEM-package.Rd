\encoding{utf8} 
\name{funFEM-package}
\alias{funFEM-package}
\docType{package}
\title{
Model-based clustering in the discriminative functional subspaces with the funFEM algorithm
}
\description{
The package provides the funFEM algorithm (Bouveyron et al., 2014) which allows to cluster functional data by modeling the curves within a common and discriminative functional subspace. 
}
\details{
\tabular{ll}{
Package: \tab funFEM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-06\cr
License: \tab GPL-2\cr
}
}
\author{
Charles Bouveyron

Maintainer: <charles.bouveyron@parisdescartes.fr>
}
\references{
C. Bouveyron, E. Côme and J. Jacques, The discriminative functional mixture model for the analysis of bike sharing systems, Preprint HAL n.01024186, University Paris Descartes, 2014.
}
\keyword{ package }
\examples{
# Clustering the well-known "Canadian temperature" data (Ramsay & Silverman)
basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4)
fdobj <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],basis,
        fdnames=list("Day", "Station", "Deg C"))$fd
res = funFEM(fdobj,K=4)

# Visualization of the partition and the group means
par(mfrow=c(1,2))
plot(fdobj,col=res$cls,lwd=2,lty=1)
fdmeans = fdobj; fdmeans$coefs = t(res$prms$my)
plot(fdmeans,col=1:max(res$cls),lwd=2)
}
