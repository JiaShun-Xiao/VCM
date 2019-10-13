LMM_gibbs <- function(X,y,Z=NULL,burnin=10,maxIter=1500,tol=1e-6,a0=0,b0=0,c0=0,d0=0,verbose=T,dispEvery=20){
  p <- ncol(X)
  n <- length(y)

  # Project covariates
  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }
  ZZ <- t(Z)%*%Z
  SZy <- solve(ZZ,t(Z)%*%y)  # This is the fixed effect beta0
  y <- y - Z %*% SZy  # project covariates, mean of y becomes 0
  SZX <- solve(ZZ,t(Z)%*%X)
  X <- X - Z %*% SZX

  Xsd <- sqrt(colMeans(X^2))
  X <- t(t(X)/Xsd)/sqrt(p)  # scale X

  xx <- colSums(X^2)
  xty <- t(y)%*%X

  # record gibbs results
  beta_all <- matrix(0,p,maxIter)
  s2_all <- matrix(0,p,maxIter)
  sb2_all <- rep(0,maxIter)
  se2_all <- rep(0,maxIter)
  PVE_all <- rep(0,maxIter)

  # Initialize
  yt <- rep(0,n)
  beta <- rep(0,p)
  se2 <- var(y)/2
  sb2 <- se2

  for(iter in 1:(burnin+maxIter)){
    # Sampling beta
    se2_sb2 <- se2/sb2
    xxsig <- xx + c(se2_sb2)
    s2 <- c(se2)/xxsig
    for(j in 1:p){
      yt_j <- yt - beta[j]*X[,j]
      mu_j <- (xty[j]-t(yt_j)%*%X[,j]) / xxsig[j]
      beta[j] <- mu_j + sqrt(s2[j])*rnorm(1)
      yt <- yt_j + beta[j]*X[,j]
    }

    # Sampling sb2
    a <- a0 + p/2
    b <- b0 + sum(beta^2)/2
    inv_sb2 <- rgamma(1,shape=a,rate=b)
    sb2 <- 1/inv_sb2

    # Sampling se2
    c <- c0 + n/2
    d <- d0 + sum((y-yt)^2)/2
    inv_se2 <- rgamma(1,shape = c,rate = d)
    se2 <- 1/inv_se2

    PVE <- var(yt)/var(y)

    # Record sampling results
    if(iter > burnin){
      beta_all[,iter-burnin] <- beta
      s2_all[,iter-burnin] <- s2
      sb2_all[iter-burnin] <- sb2
      se2_all[iter-burnin] <- se2
      PVE_all[iter-burnin] <- PVE
    }

    # Display
    if(verbose && (iter%%dispEvery)==0){
      if(iter<=burnin){
        cat("Burn-in iteration ",iter,"\n",sep="")
      } else{
        cat("Iteration ",iter-burnin,"\n",sep="")
      }
    }
  }

  mu <- rowMeans(beta_all) / Xsd / sqrt(p)
  cov_mu <- cov(t(beta_all))
  sb2 <- mean(sb2_all)
  var_sb2 <- var(sb2_all)
  se2 <- mean(se2_all)
  var_se2 <- var(se2_all)
  PVE <- mean(PVE_all)
  var_PVE <- var(PVE_all)

  H <- matrix(0,2,3)
  colnames(H) <- c("sb2","se2","PVE")
  H[1,] <- c(sb2,se2,PVE)
  H[2,] <- sqrt(c(var_sb2,var_se2,var_PVE))
  sampler <- list(beta=beta_all,s2=s2_all,sb2=sb2_all,se2=se2_all,PVE=PVE_all)

  return(list(beta0=SZy,mu=mu,cov_mu=cov_mu,H=H,sampler=sampler))
}

