get_mu <- function(X,y,Z,se2,sb2){
  p <- ncol(X)
  n <- length(y)

  ym <- mean(y)
  y  <- y-ym  # centerize y
  Xm <- colMeans(X)
  X <- scale(X,center = T,scale = F)  # centerize X
  Xsd <- sqrt(colMeans(X^2))
  X <- t(t(X)/Xsd)/sqrt(p)  # scale X
  K <- X %*% t(X)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  eigenK <- eigen(K)
  eVal <- eigenK$values
  eVec <- eigenK$vectors

  yt <- t(eVec)%*%y
  Zt <- t(eVec)%*%Z

  D <- eVal/se2 + 1/sb2
  d <- 1/(D*se2*sb2) # <=> d <- 1/(sb2*eVal + se2) in the original paper, d = 1/D/se2/sb2
  beta0 <- solve((t(Zt*d))%*%Zt) %*% (t(Zt) %*% (yt*d))

  mu <- 1/se2 * t(X) %*% (eVec %*% ((yt - Zt%*%beta0) / D))  # after covnergence obtain mu using one E-step
  # recover beta0 and mu
  mu <- mu/Xsd/sqrt(p)
  if(ncol(Z)==1){
    beta0 <- beta0 + ym - colSums(mu*Xm)
  } else{
    beta0 <- beta0 + solve(ZZ,colSums(Z)*(ym-sum(mu*Xm)))  # = beta0 - (Z^TZ)^-1Z^T[(y^bar-X^bar%*%mu),...,(y^bar-X^bar%*%mu)]
  }
}
