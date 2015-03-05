\encoding{utf8} 
\name{velov}
\alias{velov}
\docType{data}
\title{
The Vélov data set
}
\description{
This data set contains data from the bike sharing system of Lyon, called Vélo'v. The data are loading profiles of the bike stations over one week. The data were collected every hour during the period Sunday 9th March - Sunday 16th March, 2014.
}
\usage{data(velov)}
\format{
The format is:

  - data: the loading profiles (nb of available bikes / nb of bike docks) of the 345 stations at 181 times.
  
  - position: the longitude and latitude of the 345 bike stations.
  
  - dates: the download dates.
  
  - bonus: indicates if the station is on a hill (bonus = 1).
  
  - names: the names of the stations.
}
\source{
The real time data are available at https://developer.jcdecaux.com/ (with an api key).
}
\references{
The data were first used in C. Bouveyron, E. Côme and J. Jacques, The discriminative functional mixture model for the analysis of bike sharing systems, Preprint HAL n.01024186, University Paris Descartes, 2014.
}
\examples{
data(velov)
matplot(t(velov$data[1:5,]),type='l',lty=1,col=2:5,xaxt='n',lwd=2,ylim=c(0,1))
axis(1,at=seq(5,181,6),labels=velov$dates[seq(5,181,6)],las=2)
}
\keyword{datasets}
