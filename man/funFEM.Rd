\encoding{utf8} 
\name{funFEM}
\alias{funFEM}
\title{
The funFEM algorithm for the clustering of functional data.
}
\description{
The funFEM algorithm allows to cluster time series or, more generally, functional data. It is based on a discriminative functional mixture model which allows the clustering of the data in a unique and discriminative functional subspace. This model presents the advantage to be parsimonious and can therefore handle long time series.
}
\usage{
funFEM(fd, K=2:6, model = "AkjBk", crit = "bic", init = "kmeans", Tinit = c(), maxit = 50,
  eps = 1e-06, disp = FALSE, lambda = 0, graph = FALSE)
}
\arguments{
  \item{fd}{
a functional data object produced by the fda package.
}
  \item{K}{
an integer vector specifying the numbers of mixture components (clusters) among which the model selection criterion will choose the most appropriate number of groups. Default is 2:6.
}
  \item{model}{
a vector of discriminative latent mixture (DLM) models to fit. There are 12 different models: "DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB".  The option "all" executes the funFEM algorithm on the 12 models and select the best model according to the maximum value obtained by model selection criterion.
}
  \item{crit}{
the criterion to be used for model selection ('bic', 'aic' or 'icl'). 'bic' is the default.
}
  \item{init}{
the initialization type ('random', 'kmeans' of 'hclust'). 'kmeans' is the default.
}
  \item{Tinit}{
a n x K matrix which contains posterior probabilities for initializing the algorithm (each line corresponds to an individual).
}
  \item{maxit}{
the maximum number of iterations before the stop of the Fisher-EM algorithm. 
}
  \item{eps}{
the threshold value for the likelihood differences to stop the Fisher-EM algorithm.
}
  \item{disp}{
if true, some messages are printed during the clustering. Default is false.
}
  \item{lambda}{
the l0 penalty (between 0 and 1) for the sparse version. See (Bouveyron et al., 2014) for details. Default is 0.
}
  \item{graph}{
if true, it plots the evolution of the log-likelhood. Default is false.
}
}
\value{
A list is returned: 
  \item{model}{the model name.}
  \item{K}{the number of groups.}
  \item{cls}{the group membership of each individual estimated by the Fisher-EM algorithm.}
  \item{P}{the posterior probabilities of each individual for each group.}
  \item{prms}{the model parameters.}
  \item{U}{the orientation of the functional subspace according to the basis functions.}
  \item{aic}{the value of the Akaike information criterion.}
  \item{bic}{the value of the Bayesian information criterion.}
  \item{icl}{the value of the integrated completed likelihood criterion.}
  \item{loglik}{the log-likelihood values computed at each iteration of the FEM algorithm.}
  \item{ll}{the log-likelihood value obtained at the last iteration of the FEM algorithm.}
  \item{nbprm}{the number of free parameters in the model.}
  \item{call}{the call of the function.}
  \item{plot}{some information to pass to the plot.fem function.}
  \item{crit}{the model selction criterion used.}
}
\references{
C. Bouveyron, E. Côme and J. Jacques, The discriminative functional mixture model for the analysis of bike sharing systems, Preprint HAL n.01024186, University Paris Descartes, 2014.
}
\author{
Charles Bouveyron
}
\examples{
# Clustering the well-known "Canadian temperature" data (Ramsay & Silverman)
basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4)
fdobj <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],basis,
        fdnames=list("Day", "Station", "Deg C"))$fd
res = funFEM(fdobj,K=4)

# Visualization of the partition and the group means
par(mfrow=c(1,2))
plot(fdobj); lines(fdobj,col=res$cls,lwd=2,lty=1)
fdmeans = fdobj; fdmeans$coefs = t(res$prms$my)
plot(fdmeans); lines(fdmeans,col=1:max(res$cls),lwd=2)

# Visualization in the discriminative subspace (projected scores)
par(mfrow=c(1,1))
plot(t(fdobj$coefs) \%*\% res$U,col=res$cls,pch=19,main="Discriminative space")

\donttest{
###############################################################################
# Analysis of the Velib data set

# Load the velib data and smoothing
data(velib)
basis<- create.fourier.basis(c(0, 181), nbasis=25)
fdobj <- smooth.basis(1:181,t(velib$data),basis)$fd

# Clustrering with FunFEM
res = funFEM(fdobj,K=6,model='AkjBk',init='kmeans',lambda=0,disp=TRUE)

# Visualization of group means
fdmeans = fdobj; fdmeans$coefs = t(res$prms$my)
plot(fdmeans); lines(fdmeans,col=1:res$K,lwd=2,lty=1)
axis(1,at=seq(5,181,6),labels=velib$dates[seq(5,181,6)],las=2)

# # Choice of K (may be long!)
# res = funFEM(fdobj,K=2:20,model='AkjBk',init='kmeans',lambda=0,disp=TRUE)
# plot(2:20,res$plot$bic,type='b',xlab='K',main='BIC')

# Computation of the closest stations from the group means
par(mfrow=c(3,2))
for (i in 1:res$K) {
  matplot(t(velib$data[which.max(res$P[,i]),]),type='l',lty=i,col=i,xaxt='n',
          lwd=2,ylim=c(0,1))
  axis(1,at=seq(5,181,6),labels=velib$dates[seq(5,181,6)],las=2)
  title(main=paste('Cluster',i,' - ',velib$names[which.max(res$P[,i])]))
}

# Visualization in the discriminative subspace (projected scores)
par(mfrow=c(1,1))
plot(t(fdobj$coefs) \%*\% res$U,col=res$cls,pch=19,main="Discriminative space")
text(t(fdobj$coefs) \%*\% res$U)

# # Spatial visualization of the clustering (with library ggmap)
# library(ggmap)
# Mymap = get_map(location = 'Paris', zoom = 12, maptype = 'terrain')
# ggmap(Mymap) + geom_point(data=velib$position,aes(longitude,latitude),
#                           colour = I(res$cl), size = I(3))

# FunFEM clustering with sparsity
res2 = funFEM(fdobj,K=res$K,model='AkjBk',init='user',Tinit=res$P,
              lambda=0.01,disp=TRUE)

# Visualization of group means and the selected functional bases
split.screen(c(2,1))
fdmeans = fdobj; fdmeans$coefs = t(res2$prms$my)
screen(1); plot(fdmeans,col=1:res2$K,xaxt='n',lwd=2) 
axis(1,at=seq(5,181,6),labels=velib$dates[seq(5,181,6)],las=2)
basis$dropind = which(rowSums(abs(res2$U))==0)
screen(2); plot(basis,col=1,lty=1,xaxt='n',xlab='Disc. basis functions')
axis(1,at=seq(5,181,6),labels=velib$dates[seq(5,181,6)],las=2)
close.screen(all=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ functional data }% __ONLY ONE__ keyword per line
