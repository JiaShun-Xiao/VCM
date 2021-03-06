MQS <- function(file_z,file_ref,file_cov=NULL,K=NULL,X=NULL,snp_list=NULL,sd_method="LD_block",pop="EUR",file_log=""){
  if(nchar(file_log)>0) sink(file_log,append=F,split=T)
  cat("Summary statistics file: ",file_z,"\n",sep="")
  cat("Reference file: ",file_ref,"\n",sep="")
  if(! is.null(file_cov)) cat("Covariates file: ",file_cov,"\n",sep="")
  if(! is.null(file_log))cat("Writing to log file: ",file_log,"\n",sep="")

  cat("Reading data from summary statisitcs...\n")
  z_G <- read.table(file_z,header = T,stringsAsFactors = F)
  colnames(z_G) <- c("SNP","N","Z","A1","A2")
  cat(nrow(z_G)," SNPs found in summary statistics.\n",sep = "")

  cat("Reading data from reference penal...\n")
  snp_ref <- read.table(paste(file_ref,".bim",sep=""),stringsAsFactors = F)
  colnames(snp_ref) <- c("CHR","SNP","POS","BP","A1","A2")
  cat(nrow(snp_ref)," SNPs found in reference penal.\n",sep = "")

  snp_match <- intersect(snp_ref$SNP,z_G$SNP)
  if(!is.null(snp_list)) snp_match <- snp_match[snp_match%in%snp_list]

  idx_z <- match(snp_match,z_G$SNP)
  z_G <- z_G[idx_z,]

  idx_ref <- match(snp_match,snp_ref$SNP)
  snp_ref <- snp_ref[idx_ref,]
  cat(length(snp_match)," SNPs are matched in both files.\n",sep = "")

  # replace T with A, replace G with C; A=1, C=2
  z_G$A1[z_G$A1=="T"] <- "A"
  z_G$A2[z_G$A2=="T"] <- "A"
  z_G$A1[z_G$A1=="G"] <- "C"
  z_G$A2[z_G$A2=="G"] <- "C"
  z_A1 <- ifelse(z_G$A1=="A",1,2)
  z_A2 <- ifelse(z_G$A2=="A",1,2)
  # z_G$A1 <- as.factor(z_G$A1)
  # z_G$A2 <- as.factor(z_G$A2)

  snp_ref$A1[snp_ref$A1=="T"] <- "A"
  snp_ref$A2[snp_ref$A2=="T"] <- "A"
  snp_ref$A1[snp_ref$A1=="G"] <- "C"
  snp_ref$A2[snp_ref$A2=="G"] <- "C"
  snp_A1 <- ifelse(snp_ref$A1=="A",1,2)
  snp_A2 <- ifelse(snp_ref$A2=="A",1,2)
  # snp_ref$A1 <- as.factor(snp_ref$A1)
  # snp_ref$A2 <- as.factor(snp_ref$A2)

  # remove fault SNPs
  snps_rm <- (snp_A1+snp_A2)!=(z_A1+z_A2)

  snp_match <- snp_match[!snps_rm]
  idx_ref <- idx_ref[!snps_rm]
  idx_z <- idx_z[!snps_rm]
  snp_ref <- snp_ref[!snps_rm,]
  z_G <- z_G[!snps_rm,]
  cat(sum(snps_rm)," SNPs are removed because of ambiguity; ",length(snp_match)," SNPs remained.\n",sep = "")

  # calculate kinship matrix if needed
  if(is.null(K) & is.null(X)){
    cat("Calculating kinship matrix...\n")
    ref <- read_data(file_ref)
    X <- ref$X[,idx_ref]
    X[X==3] <- 0
    no_var <- apply(X,2,function(x) all(diff(x)==0))
    if(sum(no_var)!=0){
      snp_match <- snp_match[!no_var]
      X <- X[,!no_var]
      snp_ref <- snp_ref[!no_var,]
      z_G <- z_G[!no_var,]
      cat(sum(no_var)," SNPs are removed because of no variation; ",ncol(X)," SNPs remained.\n",sep = "")
    }
    X <- scale(X)
    K <- X%*%t(X)/ncol(X)
    rm(ref,X)
    gc()
  } else if(!is.null(X)){
    cat("Calculating kinship matrix...\n")
    X <- X[,idx_ref]
    X[X==3] <- 0
    no_var <- apply(X,2,function(x) all(diff(x)==0))
    if(sum(no_var)!=0){
      snp_match <- snp_match[!no_var]
      X <- X[,!no_var]
      snp_ref <- snp_ref[!no_var,]
      z_G <- z_G[!no_var,]
      cat(sum(no_var)," SNPs are removed because of no variation; ",ncol(X)," SNPs remained.\n",sep = "")
    }
    X <- scale(X)
    K <- X%*%t(X)/ncol(X)
    rm(ref,X)
    gc()
  }

  # based on the minor allel of reference penal, correct for the z score in summary statistics
  ind <- z_G$A1 != snp_ref$A1
  z_score <- z_G$Z
  z_score[ind] <- -z_score[ind]
  cat(sum(ind)," SNPs have different minor allel, z-scores are corrected according to reference penal.\n",sep = "")

  if(sd_method=="LD_block"){
    cat("Assigning SNPs to LD Blocks...\n")
    block <- read.table(system.file("extdata", paste0(pop,"_fourier_ls-all.bed"), package = "VCM"),header = T)
    group <- rep(0,nrow(z_G))
    idx_group <- 1
    for(i in 1:22){
      block_i <- block[block$chr==paste("chr",i,sep=""),]
      # chr_i <- snp_ref[snp_ref$CHR==i,]
      n_block <- nrow(block_i)

      for(j in 1:n_block){
        tmp <- with(snp_ref,CHR==i & BP>=block_i$start[j] & BP<block_i$stop[j])
        if(sum(tmp!=0)){
          group[tmp] <- idx_group
          idx_group <- idx_group+1
        }
      }
    }
  } else if(sd_method=="Jackknife"){
    group <- NULL
  } else if(sd_method=="Chromosome"){
    group <- snp_ref$CHR
  }

  cat("Calculate PVE...\n")
  if(is.null(file_cov)){
    fit <- MoM_2var_ss(K,z_score,z_G$N,Z=NULL,group=group)
  } else {
    covar <- read.table(file_cov,header = F)
    covar <- data.matrix(covar)
    fit <- MoM_2var_ss(K,z_score,z_G$N,Z=covar,group=group)
  }


  cat("Done.\n")
  cat("h2 = ",fit$h," (",fit$se_h,")\n",sep="")

  if(nchar(file_log)>0) sink()
  return(c(fit,snps=list(snp_match)))
}
