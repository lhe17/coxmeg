
#' Perform GWAS using a Cox mixed-effects model with plink files as input
#'
#' \code{coxmeg_plink} first estimates the variance component under a null model with only cov if tau is not given, and then analyzing each SNP in the plink files.
#' 
#' @section About \code{corr}:
#' The subjects in \code{corr} must be in the same order as in the plink fam file.
#' @section About missing values:
#' \code{pheno} -9 for missing values, \code{cov_file} NA for missing values.
#' @section About temporary files:
#' The function will create a temporary gds file with approximately the same size as the bed file. The temporary file will be removed when the analysis is done.
#' @section About \code{type}:
#' Specifying the type of the relatedness matrix (whether it is block-diagonal, general sparse, or dense). In the case of multiple relatedness matrices, it refers to the type of the sum of these matrices.
#' \itemize{ 
#' \item{"bd"}{ - used for a block-diagonal relatedness matrix, or a sparse matrix the inverse of which is also sparse. }
#' \item{"sparse"}{ - used for a general sparse relatedness matrix the inverse of which is not sparse.}
#' \item{"dense"}{ - used for a dense relatedness matrix.}
#' }
#' @section About \code{spd}:
#' When \code{spd=TRUE}, the relatedness matrix is treated as SPD. If the matrix is SPSD or not sure, use \code{spd=FALSE}.
#' @section About \code{solver}:
#' Specifying which method is used to solve the linear system in the optimization algorithm.  
#' \itemize{ 
#' \item{"1"}{ - Cholesky decompositon (RcppEigen:LDLT) is used to solve the linear system.}
#' \item{"2"}{ - PCG is used to solve the linear system. When \code{type='dense'}, it is recommended to set \code{solver=2} to have better computational performance.}
#' }
#' @section About \code{detap}:
#' Specifying the method to compute the log-determinant for estimating the variance component(s).
#' \itemize{ 
#' \item{"exact"}{ - the exact log-determinant is computed for estimating the variance component.}
#' \item{"diagonal"}{ - uses diagonal approximation and is only effective for a sparse relatedness matrix.}
#' \item{"slq"}{ - uses stochastic lanczos quadrature approximation. It uses the Lanczos algorithm to compute the weights and nodes. When type is 'bd' or 'sparse', it is often faster than 'gkb' and has the same accuracy. When type='dense', it is fater than 'gkb' by around half, but can be inaccurate if the relatedness matrix is (almost) singular.}
#' \item{"gkb"}{ - uses stochastic lanczos quadrature approximation. It uses the Golub-Kahan bidiagonalization algorithm to compute the weights and nodes. It is robust against an (almost) singular relatedness matrix when type='dense', but it is generally slower than 'slq'.}  
#' }
#' 
#' 
#' @param pheno A string value indicating the file name or the full path of a pheno file. The files must be in the working directory if the full path is not given. The file is in plink pheno format, containing the following four columns, family ID, individual ID, time and status. The status is a binary variable (1 for events/0 for censored).
#' @param corr A relatedness matrix or a List object of matrices if there are multiple relatedness matrices. They can be a matrix or a 'dgCMatrix' class in the Matrix package. The matrix (or the sum if there are multiple) must be symmetric positive definite or symmetric positive semidefinite. The order of subjects must be consistent with that in pheno.
#' @param type A string indicating the sparsity structure of the relatedness matrix. Should be 'bd' (block diagonal), 'sparse', or 'dense'. See details.
#' @param bed A optional string value indicating the file name or the full path of a plink bed file (without .bed). The files must be in the working directory if the full path is not given. If not provided, only the variance component will be returned.
#' @param tmp_dir A optional directory to store temporary .gds files. The directory needs to be specified when \code{bed} is provided. 
#' @param tau An optional positive value or vector for the variance component(s). If tau is given, the function will skip estimating the variance component, and use the given tau to analyze the SNPs.
#' @param cov_file An optional string value indicating the file name or the full path of a covariate file. The files must be in the working directory if the full path is not given. Same as the cov file in plink, the first two columns are family ID and individual ID. The covariates are included in the null model for estimating the variance component. The covariates can be quantitative or binary values. Categorical variables need to be converted to dummy variables.
#' @param eps An optional positive value indicating the relative convergence tolerance in the optimization algorithm. Default is 1e-6. A smaller value (e.g., 1e-8) can be used for better precision of the p-values in the situation where most SNPs under investigation have a very low minor allele count (<5).
#' @param min_tau An optional positive value indicating the lower bound in the optimization algorithm for the variance component tau. Default is 1e-4.
#' @param max_tau An optional positive value indicating the upper bound in the optimization algorithm for the variance component tau. Default is 5.
#' @param opt An optional logical scalar for the Optimization algorithm for estimating the variance component(s). Can be one of the following values: 'bobyqa', 'Brent', 'NM', or 'L-BFGS-B' (only for >1 variance components). Default is 'bobyqa'.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. See details.
#' @param detap An optional string indicating whether to use an approximation for log-determinant. Can be 'exact', 'diagonal', 'gkb', or 'slq'. Default is NULL, which lets the function select a method based on 'type' and other information. See details.
#' @param solver An optional binary value taking either 1 or 2. Default is 1. See details.
#' @param maf An optional positive value. All SNPs with MAF<maf in the bed file will not be analyzed. Default is 0.05.
#' @param score An optional logical value indicating whether to perform a score test. Default is FALSE.
#' @param threshold An optional non-negative value. If threshold>0, coxmeg_m will reestimate HRs for those SNPs with a p-value<threshold by first estimating a variant-specific variance component. Default is 0.
#' @param order An optional integer value starting from 0. Only effective when \code{dense=FALSE}. It specifies the order of approximation used in the inexact newton method. Default is NULL, which lets coxmeg choose an optimal order.
#' @param verbose An optional logical value indicating whether to print additional messages. Default is TRUE.
#' @param mc An optional integer scalar specifying the number of Monte Carlo samples used for approximating the log-determinant when \code{detap='gkb'} or \code{detap='slq'}. Default is 100.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return p: The p-value of each SNP.
#' @return tau: The estimated variance component.
#' @return rank: The rank of the relatedness matrix.
#' @return nsam: Actual sample size.
#' @keywords Cox mixed-effects model
#' @export coxmeg_plink
#' @examples
#' library(Matrix)
#' library(MASS)
#' library(coxmeg)
#' 
#' ## build a block-diagonal relatedness matrix
#' n_f <- 600
#' mat_list <- list()
#' size <- rep(5,n_f)
#' offd <- 0.5
#' for(i in 1:n_f)
#' {
#'   mat_list[[i]] <- matrix(offd,size[i],size[i])
#'   diag(mat_list[[i]]) <- 1
#' }
#' sigma <- as.matrix(bdiag(mat_list))
#' 
#' ## Estimate variance component under a null model
#' pheno = system.file("extdata", "ex_pheno.txt", package = "coxmeg")
#' cov = system.file("extdata", "ex_cov.txt", package = "coxmeg")
#' bed = system.file("extdata", "example_null.bed", package = "coxmeg")
#' bed = substr(bed,1,nchar(bed)-4)
#' re = coxmeg_plink(pheno,sigma,type='bd',bed=bed,tmp_dir=tempdir(),cov_file=cov,
#' detap='diagonal',order=1)
#' re


coxmeg_plink <- function(pheno,corr,type,bed=NULL,tmp_dir=NULL,cov_file=NULL,tau=NULL,maf=0.05,
                         min_tau=1e-04,max_tau=5,eps=1e-6,order=NULL,detap=NULL,opt='bobyqa',score=FALSE,threshold=0,
                         solver=NULL,spd=TRUE,mc=100,verbose=TRUE){
  
  if(eps<0)
  {eps <- 1e-6}
  
  slqd <- 8
  
  if(is.null(order))
  {order_t = 1}else{order_t = as.integer(order)}
  
  if(is.null(bed)==FALSE)
  {
    if((is.null(tmp_dir)==TRUE) || (file.exists(tmp_dir)==FALSE))
    {stop('The temporary directory is not specified or does not exist.')}
  }
  
  ncm <- 1
  if(is.list(corr))
  {
    ncm <- length(corr)
    if(ncm==1)
    {
      corr <- corr[[1]]
    }
  }
  
  cd = getwd()
  
  if(is.null(bed)==FALSE)
  {
    bedbimfam.fn = paste0(bed,c('.bed','.fam','.bim'))
    
    if((sum(file.exists(bedbimfam.fn))<3))
    {
      bedbimfam.fn = paste0(cd,'/',bed,c('.bed','.fam','.bim'))
      if((sum(file.exists(bedbimfam.fn))<3))
      {stop('Cannot find the genotype files.')}
    }
  }
  
  pheno.fn = paste0(cd,'/',pheno)
  pfs = c(pheno,pheno.fn)
  if(sum(file.exists(pfs))<1)
  {
    stop('Cannot find the phenotype file.')
  }else{
    phenod = read.table(pfs[file.exists(pfs)][1],header=FALSE,stringsAsFactors=FALSE) 
    nphenod <- nrow(phenod)
  }
  
  cov.fn = cov = NULL
  if(is.null(cov_file)==FALSE)
  {
    cov.fn = paste0(cd,'/',cov_file)
    pfs = c(cov_file,cov.fn)
    if(sum(file.exists(pfs))>0)
    {
      cov = read.table(pfs[file.exists(pfs)][1],header=FALSE,na.strings = "NA",stringsAsFactors=FALSE)
    }else{
      stop('Cannot find the covariate file.')
    }
  }
  
  if(is.null(cov)==TRUE)
  {samind = which((phenod[,3]!=-9)&(phenod[,4]!=-9))}else{
    covna = apply(as.matrix(cov[,3:ncol(cov)]),1,function(x) sum(is.na(x)))
    samind = which((phenod[,3]!=-9)&(phenod[,4]!=-9)&(covna==0))
  }
  
  outcome = as.matrix(phenod[samind,c(3,4)])
  samid = as.character(phenod[samind,2])
  
  if(is.null(cov)==FALSE)
  {cov = as.matrix(cov[samind,3:ncol(cov)])}
  if(ncm==1)
  {
    if(sum(dim(corr)==nphenod)<2)
    {stop('The dimension of corr is not consistent with the sample size in pheno.')}
    corr = corr[samind,samind]
  }else{
    for(i in 1:ncm)
    {
      if(sum(dim(corr[[i]])==nphenod)<2)
      {stop('The dimension of a matrix in corr is not consistent with the sample size in pheno.')}
      corr[[i]] <- corr[[i]][samind,samind,drop = FALSE]
    }
  }
  if(verbose==TRUE)
  {message(paste0(length(samind), ' subjects remain after removing subjects with missing phenotype or covariates.'))}
  
  check_input(outcome,corr,type,detap,ncm,spd)
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    samid = samid[-rem]
    outcome <- outcome[-rem,,drop = FALSE]
    if(ncm==1)
    {
      corr <- corr[-rem,-rem,drop = FALSE]
    }else{
      for(i in 1:ncm)
      {
        corr[[i]] <- corr[[i]][-rem,-rem,drop = FALSE]
      }
    }
    if(is.null(cov)==FALSE)
    {cov <- as.matrix(cov[-rem,,drop = FALSE])}
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  n <- nrow(outcome)
  
  if(ncm==1)
  {
    nz <- nnzero(corr)
  }else{
    nz <- nnzero(Reduce('+',corr))
  }
  if( nz > ((as.double(n)^2)/2) )
  {type <- 'dense'}
  
  u <- rep(0,n)
  if(is.null(cov)==FALSE)
  {
    x_sd = which(as.vector(apply(cov,2,sd))>0)
    x_ind = length(x_sd)
    if(x_ind==0)
    {
      warning("The covariates are all constants after the removal of subjects.")
      k <- 0
      beta <- numeric(0)
      cov <- matrix(0,0,0)
    }else{
      k <- ncol(cov)
      if(x_ind<k)
      {
        warning(paste0(k-x_ind," covariate(s) is/are removed because they are all constants after the removal of subjects."))
        cov = cov[,x_sd,drop=FALSE]
        k <- ncol(cov)
      }
      beta <- rep(0,k)
    }
  }else{
    k <- 0
    beta <- numeric(0)
    cov <- matrix(0,0,0)
  }
  
  d_v <- outcome[,2]
  
  ## risk set matrix
  ind <- order(outcome[,1])
  ind <- as.matrix(cbind(ind,order(ind)))
  rk <- rank(outcome[ind[,1],1],ties.method='min')
  # n1 <- sum(d_v>0)
  rs <- rs_sum(rk-1,d_v[ind[,1]])
  
  spsd = FALSE
  if(spd==FALSE)
  {spsd = TRUE}
  rk_cor = n
  if(verbose==TRUE)
  {message(paste0('There is/are ',k,' covariates. The sample size included is ',n,'.'))}
  
  eigen = TRUE
  sigma_i_s <- NULL
  if(type=='dense')
  {
    if(ncm==1)
    {
      corr <- as.matrix(corr)
      if(spsd==FALSE)
      {
        sigma_i_s = chol(corr)
        sigma_i_s = as.matrix(chol2inv(sigma_i_s))
      }else{
        sigma_i_s = eigen(corr)
        if(min(sigma_i_s$values) < -1e-10)
        {
          warning("The relatedness matrix has negative eigenvalues. Please use a positive (semi)definite matrix.")
        }
        npev <- which(sigma_i_s$values<1e-10)
        if(length(npev)>0)
        {
          sigma_i_s$values[npev] = 1e-6
        }
        sigma_i_s = sigma_i_s$vectors%*%diag(1/sigma_i_s$values)%*%t(sigma_i_s$vectors)
        # corr <- ginv(corr)
        # rk_cor = n
        # spsd = FALSE
      }
      s_d <- as.vector(diag(sigma_i_s))
    }else{
      for(i in 1:ncm)
      {corr[[i]] <- as.matrix(corr[[i]])}
    }
    
    inv <- TRUE
    
    solver <- get_solver(solver,type,verbose)
    
    detap <- set_detap_dense(detap,n,spsd,ncm)
  }else{
    if(ncm==1)
    {corr <- as(corr, 'dgCMatrix')}
    si_d = s_d = NULL
    
    if(spsd==TRUE)
    {
      minei <- rARPACK::eigs_sym(corr, 1, which = "SA")$values
      if(minei < -1e-10)
      {
        stop("The relatedness matrix has negative eigenvalues. Please use a positive (semi)definite matrix.")
      }
      if(minei<1e-10)
      {Matrix::diag(corr) <- Matrix::diag(corr) + 1e-6}
      rk_cor = n
    }
    
    if(type=='bd')
    {
      if(ncm==1)
      {
        corr <- Matrix::chol2inv(Matrix::chol(corr))
        s_d <- as.vector(Matrix::diag(corr))
        inv = TRUE
      }else{
        inv <- FALSE
      }
    }else{
      if(ncm==1)
      {s_d <- as.vector(Matrix::diag(corr))}
      inv = FALSE
    }
    
    detap <- set_detap_sparse(detap,type,verbose,ncm)
    
    solver <- get_solver(solver,type,verbose)
    
    if((inv==TRUE) & (ncm==1))
    {
      corr <- as(corr,'dgCMatrix')
    }
  }
  
  if(verbose==TRUE)
  { model_info(inv,ncm,eigen,type,solver,detap) }
  
  rad = NULL
  if(detap%in%c('slq','gkb'))
  {
    rad = rbinom(n*mc,1,0.5)
    rad[rad==0] = -1
    rad = matrix(rad,n,mc)/sqrt(n)
  }
  
  iter = NULL
  
  if(is.null(tau))
  {
    tau_e = rep(0.5,ncm) 
    if(ncm==1)
    {
      new_t = switch(
        opt,
        'bobyqa' = bobyqa(tau_e, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,det=TRUE,detap=detap,inv=inv,eps=eps,lower=min_tau,upper=max_tau,solver=solver,rad=rad,slqd=slqd),
        'Brent' = optim(tau_e, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,det=TRUE,detap=detap,inv=inv,eps=eps,lower=min_tau,upper=max_tau,method='Brent',solver=solver,rad=rad,slqd=slqd),
        'NM' = optim(tau_e, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,det=TRUE,detap=detap,inv=inv,eps=eps,method='Nelder-Mead',solver=solver,rad=rad,slqd=slqd),
        stop("The argument opt should be bobyqa, Brent or NM.")
      )
    }else{
      new_t = switch(
        opt,
        'bobyqa' = bobyqa(tau_e, mll_mm, type=type, beta=beta,u=u,s_d=s_d,corr=corr,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,detap=detap,eps=eps,lower=rep(min_tau,ncm),upper=rep(max_tau,ncm),solver=solver,rad=rad,slqd=slqd),
        'L-BFGS-B' = optim(tau_e, mll_mm, type=type, beta=beta,u=u,s_d=s_d,corr=corr,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,detap=detap,eps=eps,lower=rep(min_tau,ncm),upper=rep(max_tau,ncm),method='L-BFGS-B',solver=solver,rad=rad,slqd=slqd),
        stop("The argument opt should be bobyqa or L-BFGS-B.")
      )
    }
    
    if(opt=='bobyqa')
    {iter <- new_t$iter}else{
      iter <- new_t$counts
    }
    
    tau_e <- new_t$par
    check_tau(tau_e,min_tau,max_tau,ncm)
    tau <- tau_e
  }else{
    if(min(tau) < 0)
    {stop("The variance component must be positive.")}
    if(length(tau)!=ncm)
    {stop("The number of tau must be consistent with the number of correlation matrices.")}
    tau_e = tau
  }
  
  if(is.null(bed)==TRUE)
  {
    res = list(tau=tau_e,iter=iter,rank=rk_cor,nsam=n)
    return(res)
  }
  
  if(verbose==TRUE)
  {message(paste0('The variance component is estimated. Start analyzing SNPs...'))}
  
  gds.fn = paste0("tmp_",floor(runif(1,0,1)*1e8),".gds")
  gds.fn = file.path(tmp_dir,gds.fn)
  snpgdsBED2GDS(bedbimfam.fn[1], bedbimfam.fn[2], bedbimfam.fn[3], gds.fn,verbose=verbose)
  genofile <- snpgdsOpen(gds.fn,allow.duplicate=TRUE)
  # snp_info <- snpgdsSNPRateFreq(genofile,with.id=TRUE)
  snp_ind = snpgdsSelectSNP(genofile,sample.id=samid,maf=maf,missing.rate=0,remove.monosnp=TRUE)
  
  gdssamid <- snpgdsSummary(genofile, show=FALSE)$sample.id
  gdssamid <- gdssamid[gdssamid%in%samid]
  if(min(gdssamid==samid)==0)
  {stop('The subject IDs in the genotype files are not consistent with those in the phenotype file.')}
  
  if(ncm>1)
  {
    for(i in 1:ncm)
    {
      corr[[i]] <- corr[[i]]*tau_e[i]
    }
    corr <- Reduce('+',corr)
    
    if(type=='dense')
    {
      sigma_i_s = chol(corr)
      sigma_i_s = as.matrix(chol2inv(sigma_i_s))
      s_d <- as.vector(diag(sigma_i_s))
    }else{
      corr <- as(corr,'dgCMatrix')
      if(inv==TRUE)
      {
        corr <- Matrix::chol2inv(Matrix::chol(corr))
        corr <- as(corr,'dgCMatrix')
      }
      s_d <- as.vector(Matrix::diag(corr))
    }
    tau_e <- 1
  }
  
  nsnp = length(snp_ind)
  blocks = max(100,ceiling(5e6/n))
  nblock = ceiling(nsnp/blocks)
  remain = nsnp%%blocks
  sp = blocks*((1:nblock)-1) + 1
  ep = blocks*(1:nblock)
  if(remain>0)
  {
    ep[length(ep)] = nsnp
  }
  
  # genofile <- snpgdsOpen(gds.fn,allow.duplicate=TRUE)
  
  if(score==FALSE)
  {
    sumstats <- data.frame(index=snp_ind,beta=rep(NA,nsnp),HR=rep(NA,nsnp),sd_beta=rep(NA,nsnp),p=rep(NA,nsnp))
    beta <- rep(1e-16,k+1)
    u <- rep(0,n)
    
    if(type!='dense')
    {
      if(is.null(order))
      {
        x_test = matrix(sample(c(rep(1,floor(n/2)),rep(0,ceiling(n/2))),replace = FALSE),n,1)
        est_t = microbenchmark(irls_fast_ap(0, u, tau_e, s_d, corr, inv,x_test, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,1,det=FALSE,detap=detap,solver=solver),times=2)
        est_t = mean(est_t$time)
        for(ord_o in 2:10)
        {
          est_c = microbenchmark(irls_fast_ap(0, u, tau_e, s_d, corr, inv,x_test, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,ord_o,det=FALSE,detap=detap,solver=solver),times=2)
          est_c = mean(est_c$time)
          if(est_c<(0.9*est_t))
          {
            est_t = est_c
            order_t = ord_o
          }else{
            break
          }
        }
      }
      if(verbose==TRUE)
      {message(paste0('The order is set to be ', order_t,'.'))}
    }
    
    for(bi in 1:nblock)
    {  
      snp_t = sp[bi]:ep[bi]
      X = snpgdsGetGeno(genofile, sample.id=samid,snp.id=snp_ind[snp_t],with.id=TRUE,verbose=FALSE)
      X = X$genotype
      storage.mode(X) <- 'numeric'
      
      p <- ncol(X)
      c_ind <- c()
      if(k>0)
      {
        X <- cbind(X,cov)
        c_ind <- (p+1):(p+k)
      }
      
      if(type=='dense')
      {
        cme_re <- sapply(1:p, function(i)
        {
          tryCatch({
          res <- irls_ex(beta, u, tau_e, s_d, corr, sigma_i_s,as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=detap,solver=solver)
          c(res$beta[1],res$v11[1,1])},
          warning = function(war){
            message(paste0('The estimation may not converge for SNP ',i, ' in block ',bi))
            c(NA,NA)
            },
          error = function(err){
            message(paste0('The estimation failed for SNP ',i, ' in block ',bi))
            c(NA,NA)
            }
          )
        })
      }else{
        cme_re <- sapply(1:p, function(i)
        {
          tryCatch({
          res <- irls_fast_ap(beta, u, tau_e, s_d, corr, inv,as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order_t,det=FALSE,detap=detap,solver=solver)
          c(res$beta[1],res$v11[1,1])},
          warning = function(war){
            message(paste0('The estimation may not converge for SNP ',i, ' in block ',bi))
            c(NA,NA)
            },
          error = function(err){
            message(paste0('The estimation failed for SNP ',i, ' in block ',bi))
            c(NA,NA)
            }
          )
        })
      }
      
      sumstats$beta[snp_t] <- cme_re[1,]
      sumstats$HR[snp_t] <- exp(cme_re[1,])
      sumstats$sd_beta[snp_t] <- sqrt(cme_re[2,])
      sumstats$p[snp_t] <- pchisq(sumstats$beta[snp_t]^2/cme_re[2,],1,lower.tail=FALSE)
    }
    
    top = which(sumstats$p<threshold)
    if(length(top)>0)
    {
      if(ncm==1)
      {
        if(verbose==TRUE)
        {message(paste0('Finish analyzing SNPs. Start analyzing significant SNPs using a variant-specific variance component...'))}
        
        X = snpgdsGetGeno(genofile, sample.id=samid,snp.id=snp_ind[top],with.id=TRUE,verbose=FALSE)
        X = X$genotype
        storage.mode(X) <- 'numeric'
        p <- ncol(X)
        c_ind <- c()
        if(k>0)
        {
          X <- cbind(X,cov)
          c_ind <- (p+1):(p+k)
        }
        
        # tau = tau_e
        cme_re <- sapply(1:p, function(i)
        {
          new_t = switch(
            opt,
            'bobyqa' = bobyqa(tau, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=as.matrix(X[,c(i,c_ind)]), d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,det=TRUE,detap=detap,inv=inv,eps=eps,lower=min_tau,upper=max_tau,solver=solver,rad=rad,slqd=slqd),
            'Brent' = optim(tau, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=as.matrix(X[,c(i,c_ind)]), d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,det=TRUE,detap=detap,inv=inv,eps=eps,lower=min_tau,upper=max_tau,method='Brent',solver=solver,rad=rad,slqd=slqd),
            'NM' = optim(tau, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=as.matrix(X[,c(i,c_ind)]), d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order_t,det=TRUE,detap=detap,inv=inv,eps=eps,method='Nelder-Mead',solver=solver,rad=rad,slqd=slqd)
          )
          tau_s <- new_t$par
          if(type=='dense')
          {
            res <- irls_ex(beta, u, tau_s, s_d, corr, sigma_i_s,as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=detap,solver=solver)
          }else{
            res <- irls_fast_ap(beta, u, tau_s, s_d, corr, inv,as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order_t,det=FALSE,detap=detap,solver=solver)
          }
          c(tau_s,res$beta[1],res$v11[1,1])
        })
        
        sumstats$tau = sumstats$beta_exact = sumstats$sd_beta_exact = sumstats$p_exact = NA
        sumstats$tau[top] = cme_re[1,]
        sumstats$beta_exact[top] = cme_re[2,]
        sumstats$sd_beta_exact[top] = sqrt(cme_re[3,])
        sumstats$p_exact[top] <- pchisq(sumstats$beta_exact[top]^2/cme_re[3,],1,lower.tail=FALSE)
      }else{
        if(verbose==TRUE)
        {message('Variable-specific re-estimation of the variance components is not supported for multiple correlation matrices.')}
      }
    }
  }else{
    beta <- rep(1e-16,k)
    u <- rep(0,n)
    
    if(type=='dense')
    {
      model_n <- irls_ex(beta, u, tau_e, s_d, corr, sigma_i_s,cov, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,det=FALSE,detap=detap,solver=solver)
    }else{
      model_n <- irls_fast_ap(beta, u, tau_e, s_d, corr, inv, cov, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order_t,det=FALSE,detap=detap,solver=solver)
    }
    u <- model_n$u
    if(k>0)
    {
      gamma <- model_n$beta
      eta_v <- cov%*%gamma+u
    }else{
      eta_v <- u
    }
    
    w_v <- as.vector(exp(eta_v))
    s <- as.vector(cswei(w_v,rs$rs_rs-1,ind-1,1))
    a_v <- d_v/s
    a_v_p <- a_v[ind[,1]]
    a_v_2 <- as.vector(a_v_p*a_v_p)
    a_v_p <- a_v_p[a_v_p>0]
    b_v <- as.vector(cswei(a_v,rs$rs_cs-1,ind-1,0))
    bw_v <- as.vector(w_v)*as.vector(b_v)
    deriv <- d_v - bw_v
    
    dim_v = k + n
    brc = (k+1):dim_v
    v = matrix(NA,dim_v,dim_v)
    v[brc,brc] = -wma_cp(w_v,rs$rs_cs_p-1,ind-1,a_v_p)
    diag(v[brc,brc]) = diag(v[brc,brc]) + bw_v
    if(k>0)
    {
      temp = t(cov)%*%v[brc,brc]
      v[1:k,1:k] = temp%*%cov
      v[1:k,brc] = temp
      v[brc,1:k] = t(temp)
    }
    #v[brc,brc] = v[brc,brc] + as.matrix(sigma_i_s/tau_e)
    if(type=='dense')
    {
      v[brc,brc] = v[brc,brc] + as.matrix(sigma_i_s/tau_e)
    }else{
      if(inv==FALSE)
      {corr <- Matrix::chol2inv(Matrix::chol(corr))}
      v[brc,brc] = v[brc,brc] + as.matrix(corr/tau_e)
    }
    v = chol(v)
    v = chol2inv(v)
    
    sumstats <- data.frame(index=snp_ind,score=rep(NA,nsnp),score_test=rep(NA,nsnp),p=rep(NA,nsnp))
    
    for(bi in 1:nblock)
    {  
      snp_t = sp[bi]:ep[bi]
      X = snpgdsGetGeno(genofile, sample.id=samid,snp.id=snp_ind[snp_t],with.id=TRUE,verbose=FALSE)
      X = X$genotype
      storage.mode(X) <- 'numeric'
      
      t_st <- score_test(deriv,bw_v,w_v,rs$rs_rs-1,rs$rs_cs-1,rs$rs_cs_p-1,ind-1,a_v_p,a_v_2,tau_e,v,cov,X)
      pv <- pchisq(t_st[,2],1,lower.tail=FALSE)
      sumstats$score[snp_t]=t_st[,1]
      sumstats$score_test[snp_t]=t_st[,2]
      sumstats$p[snp_t]=pv
    }
  }
  
  snplist = snpgdsSNPList(genofile)
  af_inc = snpgdsSNPList(genofile, sample.id=samid)
  snplist = cbind(snplist,af_inc[,'afreq'])
  colnames(snplist)[ncol(snplist)] = 'afreq_inc'
  snplist = snplist[match(snp_ind,snplist[,1]),]
  sumstats = cbind(snplist,sumstats)
  snpgdsClose(genofile)
  file.remove(gds.fn)
  res = list(summary=sumstats,tau=tau,rank=rk_cor,nsam=n)
  return(res)
}

