.criteria <-
function(loglik,T,prms,n){
  K = prms$K
  p = prms$p
  comp = switch(prms$model,
                'DkBk' = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K^2*(K-1)/2 + K,
                'DkB'  = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K^2*(K-1)/2 + 1,
                'DBk'  = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K*(K-1)/2 + K,
                'DB'   = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K*(K-1)/2 + 1,
                'AkjBk'= (K-1) + K*(K-1) + (K-1)*(p-K/2) + K^2,
                'AkjB' = (K-1) + K*(K-1) + (K-1)*(p-K/2) + K*(K-1)+1,
                'AkBk' = (K-1) + K*(K-1) + (K-1)*(p-K/2) + 2*K,
                'AkB'  = (K-1) + K*(K-1) + (K-1)*(p-K/2) + K+1,
                'AjBk' = (K-1) + K*(K-1) + (K-1)*(p-K/2) + (K-1)+K,
                'AjB'  = (K-1) + K*(K-1) + (K-1)*(p-K/2) + (K-1)+1,
                'ABk'  = (K-1) + K*(K-1) + (K-1)*(p-K/2) + K+1,
                'AB'   = (K-1) + K*(K-1) + (K-1)*(p-K/2) + 2)
  aic = loglik - comp # AIC criterion
  bic = loglik - 1/2 * comp * log(n) # BIC criterion
  T[T<1e-6] = 1e-6
  icl = loglik - 1/2 *  comp * log(n) - sum(T*log(T)) # ICL criterion
  list(aic=aic,bic=bic,icl=icl,nbprm=comp)
}
.estep <-
function(prms,fd,U){
  # Initialization
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  K = prms$K
  mu = prms$mean
  prop = prms$prop
  D = prms$D
  d = K-1
  
  QQ = matrix(NA,n,K)
  QQ2 = matrix(NA,n,K)
  T = matrix(NA,n,K)
  
  # Compute posterior probabilities
  for (k in 1:K){
    bk = D[k,p,p]
    mY = prms$my[k,]
    YY = Y-t(matrix(rep(mY,n),p,n)) 
    projYY = YY %*% U %*% t(U) # proj dans l'espace latent
    
    if (d==1){
      for (i in 1:n){QQ[i,k] =  1/D[k,1,1] * sum(projYY[i,]^2) + 1/D[k,p,p]*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)}
    }
    else{
      sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
      for (i in 1:n){QQ[i,k] =  projYY[i,] %*% sY %*% as.matrix(projYY[i, ],p,1) + 1/bk*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)}
    }
  }
  A = -1/2 * QQ
  loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
  for (k in 1:K) {
    # browser()
    T[,k] = 1 / rowSums(exp(0.5*(QQ[,k]*matrix(1,n,K)-QQ)))}
  list(T=T,loglik=loglik)
}
.fstep <-
function(fd,T,lambda){
  if (min(colSums(T)) <= 1) stop("One cluster is almost empty!") 
  G = t(fd$coefs)
  n = nrow(G)
  p = ncol(G)
  d = ncol(T) - 1
  basisobj <- fd$basis
  #W <- (basisobj + t(basisobj))/2
  W <- inprod(basisobj,basisobj)
  Ttilde  = t(apply(T,1,'/',sqrt(colSums(T))))
  eig <- svd(ginv(t(G)%*%G%*%W) %*% (t(G)%*%Ttilde%*%t(Ttilde)%*%G%*%W),nu=d,nv=0)
  U   = eig$u[,1:d]

  # Sparse version
  if (lambda>0){
    X = G %*% U
    Utilde = U
    for (i in 1:d){ x.predict = X[,i]
        res.enet = enet(G,x.predict,intercept=FALSE)
        coef     = predict.enet(res.enet,G,type="coefficients",mode="fraction",s=lambda)$coef
        Utilde[,i] = coef/ sqrt(sum(coef^2))
    }
    U = svd(Utilde)$u
  }
  U
}
.FunFEM.main <-
function(fd,K,model='AkjBk',init='kmeans',lambda=0,Tinit=c(),maxit=50,eps=1e-8,graph=F){
  # Initialization
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  # 
  # New objects
  Lobs = rep(c(-Inf),1,(maxit+1))
  # 	
  # Initialization of T
  if (init=='user'){ T = Tinit}
  else if (init=='kmeans'){
    T = matrix(0,n,K)
    ind = kmeans(Y,K)$cluster
    for (k in 1:K){ T[which(ind==k),k] = 1}
  }
  else if (init=='random'){
    T = t(rmultinom(n,1,c(rep(1/K,K))))
  }
  else if (init=='hclust'){
    T   = matrix(0,n,K)
    ind = cutree(hclust(dist(Y),method='ward.D2'),K)
    for (k in 1:K){ T[which(ind==k),k] = 1}
  }
  
  V         = .fstep(fd,T,lambda)
  prms      = .mstep(fd,V,T,model=model)
  res.estep = .estep(prms,fd,V)
  T         = res.estep$T
  Lobs[1]   = res.estep$loglik
  
  # Main loop
  Linf_new  = Lobs[1]
  for (i in 1:maxit){
    # The three main steps F, M and E
    #cat('.')
    V         = .fstep(fd,T,lambda)
    prms      = .mstep(fd,V,T,model=model)
    res.estep = .estep(prms,fd,V)
    T         = res.estep$T
    Lobs[i+1] = res.estep$loglik
    
    # Stop criterion
    if (i>=2){
      acc      = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new = Lobs[i] + 1/(1-acc) * (Lobs[i+1] - Lobs[i])
      if (abs(Linf_new - Linf_old) < eps | is.na(Linf_new)) {break}
    }
  }
  
  # Graphic option
  if (graph){
    par(mfrow=c(1,2))
    plot(as.data.frame(as.matrix(Y) %*% V[,1:2]),col=max.col(T),xlab='axis 1',ylab='axis 2',pch=20)
    plot(Lobs[1:i],xlab='iterations',ylab='vraisemblance Espace observe',col=2,pch=20)
  }
  
  # Returning the results
  cls  = max.col(T)
  crit = .criteria(Lobs[(i+1)],T,prms,n)
  W = inprod(fd$basis,fd$basis)
  U = t(W) %*% V
  res  = list(model=model,K=K,cls=cls,P=T,prms=prms,U=U,aic=crit$aic,bic=crit$bic,icl=crit$icl,
              loglik=Lobs[2:(i+1)],ll=Lobs[i+1],nbprm=crit$nbprm)
}
.mstep <-
function(fd,U,T,model){
  # 12 different submodels: [DkBk] ... [AkjBk]
  # Initialization
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  K = ncol(T)
  d = K-1
  
  mu   = matrix(NA,K,K-1)
  m   = matrix(NA,K,p)
  prop = rep(c(NA),1,K)
  D = array(0,c(K,p,p))
  
  # Projection
  W = inprod(fd$basis,fd$basis)
  U = t(W) %*% U
  X = Y %*% U
  
  # Estimation
  for (k in 1:K){
    
    nk  = sum(T[,k])
    # Prior Probability
    prop[k] = nk / n
    # Mean in the latent space
    mu[k,]  = colSums((T[,k]*matrix(1,n,d))* X) / nk
    # Observed space
    m[k,]  = colSums(T[,k]*matrix(1,n,p)* Y) / nk
    YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
    Ck  = crossprod(T[,k]*matrix(1,n,p)* YY, YY) / (nk-1) #crossprod(x,y) = t(x) %*% y
    C   = cov(Y)
    
    # Estimation of Delta k amongst 8 submodels
    if (model=='DkBk'){
      D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='DkB'){
      D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='DBk'){
      D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='DB'){
      D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkjBk'){
      if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
        D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkjB'){
      if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
        D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkBk'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkB'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AjBk'){
      if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else {
        D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='AjB'){
      if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else{
        D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='ABk'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='AB'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))	
    }
  }
  prms = list(K=K,p=p,mean=mu,my=m,prop=prop,D=D,model=model)
}