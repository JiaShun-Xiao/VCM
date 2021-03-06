corr_MoM <- function(y,z,X1,X2,Z1=NULL,Z2=NULL){
  n1 <- length(y)
  p1 <- ncol(X1)
  n2 <- length(z)
  p2 <- ncol(X2)


  X1m <- colMeans(X1)
  ym <- mean(y)
  X1 <- scale(X1,center = T,scale = F)  # centerize X
  X1sd <- sqrt(colMeans(X1^2))
  X1 <- t(t(X1)/X1sd)/sqrt(p1)  # scale X
  K1 <- X1 %*% t(X1)

  X2m <- colMeans(X2)
  zm <- mean(z)
  X2 <- scale(X2,center = T,scale = F)  # centerize X
  X2sd <- sqrt(colMeans(X2^2))
  X2 <- t(t(X2)/X2sd)/sqrt(p2)  # scale X
  K2 <- X2 %*% t(X2)

  K12 <- X1%*%t(X2)

  y0 <- c(y,z)

  if(is.null(Z1)){
    Z1 <- matrix(1,n1,1)
    M1 <- diag(n1) - matrix(1/n1,n1,n1)
    y <- y - mean(y)
    MK1 <- K1
    MK12 <- K12
  } else{
    Z1 <- cbind(1,Z1)
    M1 <- diag(n1) - Z1%*%solve(t(Z1)%*%Z1)%*%t(Z1)
    y <- M1 %*% y
    MK1 <- M1%*%K1
    MK12 <- M1%*%K12
  }
  q1 <- ncol(Z1)

  if(is.null(Z2)){
    Z2 <- matrix(1,n2,1)
    M2 <- diag(n2) - matrix(1/n2,n2,n2)
    z <- z - mean(z)
    MK2 <- K2
    K12M <- K12
  } else{
    Z2 <- cbind(1,Z2)
    M2 <- diag(n2) - Z2%*%solve(t(Z2)%*%Z2)%*%t(Z2)
    z <- M2 %*% z
    MK2 <- M2%*%K2
    K12M <- K12%*%M2
  }
  q2 <- ncol(Z2)

  trK1 <- sum(diag(MK1))
  trK2 <- sum(diag(MK2))

  trK1sq <- sum(MK1^2)
  trK2sq <- sum(MK2^2)
  trK12sq <- sum(diag(MK12%*%t(K12M)))

  # Xy <- t(X1)%*%y
  # yKy <- sum((Xy)^2)
  # Xz <- t(X2)%*%z
  # zKz <- sum((Xz)^2)
  #
  # yKz <- sum(Xy*Xz)

  # solve for expression level vc model
  Sy <- matrix(c(trK1sq, trK1, trK1, n1-q1),2,2)
  cy <- c(t(y)%*%K1%*%y, sum(y^2))

  invSy <- solve(Sy)
  sigmay <- invSy %*% cy

  # solve for phenotype level vc model
  Sz <- matrix(c(trK2sq, trK2, trK2, n2-q2),2,2)
  cz <- c(t(z)%*%K2%*%z, sum(z^2))

  invSz <- solve(Sz)
  sigmaz <- invSz %*% cz

  # estimate covariance
  cov <- drop(t(y)%*%K12%*%z / trK12sq)

  # estimate correlation coefficient
  rho <- cov/sqrt(sigmay[1]*sigmaz[1])

  # calculate fixed effect
  Omega <- matrix(0,n1+n2,n1+n2)
  Omega[1:n1,1:n1] <- sigmay[1]*K1
  diag(Omega[1:n1,1:n1]) <- diag(Omega[1:n1,1:n1]) + sigmay[2]
  Omega[1:n1,(n1+1):(n1+n2)] <- cov * K12
  Omega[(n1+1):(n1+n2),1:n1] <- cov * t(K12)
  Omega[(n1+1):(n1+n2),(n1+1):(n1+n2)] <- sigmaz[1]*K2
  diag(Omega[(n1+1):(n1+n2),(n1+1):(n1+n2)]) <- diag(Omega[(n1+1):(n1+n2),(n1+1):(n1+n2)]) + sigmaz[2]

  inv_Omega <- chol2inv(chol(Omega))

  Z <- matrix(0,n1+n2,q1+q2)
  Z[1:n1,1:q1] <- Z1
  Z[(n1+1):(n1+n2),(q1+1):(q1+q2)] <- Z2
  # recover beta0 and mu
  beta0 <- solve(t(Z)%*%inv_Omega%*%Z,t(Z)%*%(inv_Omega%*%y0))
  mu <- t(rbind(sigmay[1]*X1,cov*X2))%*%(inv_Omega%*%(y0-Z%*%beta0))

  # re-scale beta0 and mu
  mu <- mu/X1sd/sqrt(p1)
  beta0[1] <- beta0[1] - colSums(mu*X1m)

  # Sandwich estimator for variance of sig_y and sig_1
  covBy <- matrix(0,2,2)
  if(q1==1){
    Sigmay <- sigmay[1]*K1 + sigmay[2]*M1
  } else{
    Sigmay <- sigmay[1]*MK1%*%M1 + sigmay[2]*M1
  }
  K1S <- MK1%*%Sigmay

  covBy[1,1] <- sum(K1S^2) * 2
  covBy[2,2] <- sum(Sigmay^2) * 2
  covBy[1,2] <- covBy[2,1] <- sum(K1S*Sigmay) * 2
  covSigy <- invSy %*% covBy %*% invSy

  # Sandwich estimator for variance of sig_z and sig_2
  covBz <- matrix(0,2,2)
  if(q2==1){
    Sigmaz <- sigmaz[1]*K2 + sigmaz[2]*M2
  } else{
    Sigmaz <- sigmaz[1]*MK2%*%M2 + sigmaz[2]*M2
  }
  K2S <- MK2%*%Sigmaz

  covBz[1,1] <- sum(K2S^2) * 2
  covBz[2,2] <- sum(Sigmaz^2) * 2
  covBz[1,2] <- covBz[2,1] <- sum(K2S*Sigmaz) * 2
  covSigz <- invSz %*% covBz %*% invSz

  # sandwich estimator for covariance cov
  Omegat <- Xt <- matrix(0,n1+n2,n1+n2)
  Omegat[1:n1,1:n1] <- Sigmay
  Omegat[1:n1,(n1+1):(n1+n2)] <- cov * MK12
  Omegat[(n1+1):(n1+n2),1:n1] <- cov * t(MK12)
  Omegat[(n1+1):(n1+n2),(n1+1):(n1+n2)] <- Sigmaz
  Xt[1:n1,(n1+1):(n1+n2)] <- MK12
  Xt[(n1+1):(n1+n2),1:n1] <- t(MK12)
  OmegaXt <- Omegat %*% Xt

  varByz <- sum(OmegaXt^2) / 2
  var_cov <- varByz/trK12sq/trK12sq

  # calculate variance of rho using Delta method
  g <- -c(1/sqrt(sigmay[1]*sigmaz[1]),0.5*cov/sqrt(sigmay[1]^3*sigmaz[1]),0.5*cov/sqrt(sigmay[1]*sigmaz[1]^3))
  var_rho <- sum(g^2*c(var_cov,covSigy[1,1],covSigz[1,1]))

  Sigma <- matrix(c(sigmay,sigmaz,sqrt(diag(covSigy)),sqrt(diag(covSigz))),2,4,byrow = T)

  ret <- list(Sigma=Sigma, rho=c(rho,sqrt(var_rho)), K1=K1, K2=K2, K12, covSigy=covSigy, covSigz=covSigz, cov = c(cov,sqrt(var_cov)),beta0=beta0[1:q1],mu=mu)

  class(ret) <- c("VCMgc","MoM")
  ret
}
