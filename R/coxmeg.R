
#' Fit a Cox mixed-effects model
#'
#' \code{coxmeg} returns the estimates of the variance component(s) along with the HRs and p-values of the predictors.
#' 
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
#' \item{"slq"}{ - uses stochastic lanczos quadrature approximation. It uses the Lanczos algorithm to compute the weights and nodes. When type is 'bd' or 'sparse', it is often faster than 'gkb' and has the same accuracy. When type='dense', it is faster than 'gkb' by around half, but can be inaccurate if the relatedness matrix is (almost) singular.}
#' \item{"gkb"}{ - uses stochastic lanczos quadrature approximation. It uses the Golub-Kahan bidiagonalization algorithm to compute the weights and nodes. It is robust against an (almost) singular relatedness matrix when type='dense', but it is generally slower than 'slq'.}  
#' }
#' 
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for events / 0 for censored).
#' @param corr A relatedness matrix or a List object of matrices if there are multiple relatedness matrices. They can be a matrix or a 'dgCMatrix' class in the Matrix package. The matrix (or the sum if there are multiple) must be symmetric positive definite or symmetric positive semidefinite. The order of subjects must be consistent with that in outcome.
#' @param type A string indicating the sparsity structure of the relatedness matrix. Should be 'bd' (block diagonal), 'sparse', or 'dense'. See details.
#' @param X An optional matrix of the preidctors with fixed effects. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns. 
#' @param eps An optional positive scalar indicating the relative convergence tolerance in the optimization algorithm. Default is 1e-6.
#' @param min_tau An optional positive scalar indicating the lower bound in the optimization algorithm for the variance component \code{tau}. Default is 1e-4.
#' @param max_tau An optional positive scalar indicating the upper bound in the optimization algorithm for the variance component \code{tau} Default is 5.
#' @param opt An optional logical scalar for the Optimization algorithm for estimating the variance component(s). Can be one of the following values: 'bobyqa', 'Brent', 'NM', or 'L-BFGS-B' (only for >1 variance components). Default is 'bobyqa'.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. See details.
#' @param detap An optional string indicating whether to use an approximation for log-determinant. Can be 'exact', 'diagonal', 'gkb', or 'slq'. Default is NULL, which lets the function select a method based on 'type' and other information. See details.
#' @param solver An optional bianry value that can be either 1 (Cholesky Decomposition using RcppEigen) or 2 (PCG). Default is NULL, which lets the function select a solver. See details.
#' @param order An optional integer starting from 0. Only valid when \code{dense=FALSE}. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @param verbose An optional logical scalar indicating whether to print additional messages. Default is TRUE.
#' @param mc An optional integer scalar specifying the number of Monte Carlo samples used for approximating the log-determinant when \code{detap='gkb'} or \code{detap='slq'}. Default is 100.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return p: The p-value.
#' @return iter: The number of iterations until convergence.
#' @return tau: The estimated variance component.
#' @return int_ll: The marginal likelihood (-2*log(lik)) of tau evaluated at the estimate of tau.
#' @return rank: The rank of the relatedness matrix.
#' @return nsam: Actual sample size.
#' @keywords Cox mixed-effects model
#' @export coxmeg
#' @examples
#' library(Matrix)
#' library(MASS)
#' library(coxmeg)
#' 
#' ## simulate a block-diagonal relatedness matrix
#' tau_var <- 0.2
#' n_f <- 100
#' mat_list <- list()
#' size <- rep(10,n_f)
#' offd <- 0.5
#' for(i in 1:n_f)
#' {
#'   mat_list[[i]] <- matrix(offd,size[i],size[i])
#'   diag(mat_list[[i]]) <- 1
#' }
#' sigma <- as.matrix(bdiag(mat_list))
#' n <- nrow(sigma)
#' 
#' ## simulate random effects and outcomes
#' x <- mvrnorm(1, rep(0,n), tau_var*sigma)
#' myrates <- exp(x-1)
#' y <- rexp(n, rate = myrates)
#' cen <- rexp(n, rate = 0.02 )
#' ycen <- pmin(y, cen)
#' outcome <- cbind(ycen,as.numeric(y <= cen))
#' 
#' ## fit a Cox mixed-effects model
#' re = coxmeg(outcome,sigma,type='bd',detap='diagonal',order=1)
#' re

coxmeg <- function(outcome,corr,type,X=NULL,eps=1e-6, min_tau=1e-04,max_tau=5,order=1,detap=NULL,opt='bobyqa',solver=NULL,spd=TRUE,verbose=TRUE, mc=100){
  
  if(eps<0)
  {eps <- 1e-6}
  
  slqd <- 8
  
  if(is.null(X)==FALSE)
  {
    X <- as.matrix(X)
    storage.mode(X) <- 'numeric'
  }
  outcome <- as.matrix(outcome)
  
  ncm <- 1
  if(is.list(corr))
  {
    ncm <- length(corr)
    if(ncm==1)
    {
      corr <- corr[[1]]
    }
  }
  tau <- rep(0.5,ncm)
  
  check_input(outcome,corr,type,detap,ncm,spd)
  
  nro = nrow(outcome)
  if(is.null(X)==FALSE)
  {
    if(nro!=nrow(X))
    {stop("The phenotype and predictor matrices have different sample sizes.")}
  }
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem,,drop = FALSE]
    if(is.null(X)==FALSE)
    {
      X <- as.matrix(X[-rem,,drop = FALSE])
    }
    if(ncm==1)
    {
      corr <- corr[-rem,-rem,drop = FALSE]
    }else{
      for(i in 1:ncm)
      {
        corr[[i]] <- corr[[i]][-rem,-rem,drop = FALSE]
      }
    }
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
  if(is.null(X)==FALSE)
  {
    x_sd = which(as.vector(apply(X,2,sd))>0)
    x_ind = length(x_sd)
    if(x_ind==0)
    {
      warning("The predictors are all constants after the removal of subjects.")
      k <- 0
      beta <- numeric(0)
      X <- matrix(0,0,0)
    }else{
      k <- ncol(X)
      if(x_ind<k)
      {
        warning(paste0(k-x_ind," predictor(s) is/are removed because they are all constants after the removal of subjects."))
        X = X[,x_sd,drop=FALSE]
        k <- ncol(X)
      }
      beta <- rep(0,k)
    }
  }else{
    k <- 0
    beta <- numeric(0)
    X <- matrix(0,0,0)
  }
  
  d_v <- outcome[,2]
  
  ## risk set matrix
  ind <- order(outcome[,1])
  ind <- as.matrix(cbind(ind,order(ind)))
  mode(ind) <- "integer"
  rk <- as.integer(rank(outcome[ind[,1],1],ties.method='min') - 1)
  # n1 <- sum(d_v>0)
  rs <- rs_sum(rk,d_v[ind[,1]])
  
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
        
      }
      
      s_d <- as.vector(diag(sigma_i_s))
    }else{
      for(i in 1:ncm)
      {corr[[i]] = as.matrix(corr[[i]])}
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
  
  marg_ll = 0
  if(ncm==1)
  {
    new_t = switch(
      opt,
      'bobyqa' = bobyqa(tau, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,eps=eps,lower=min_tau,upper=max_tau,solver=solver,rad=rad,slqd=slqd),
      'Brent' = optim(tau, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,eps=eps,lower=min_tau,upper=max_tau,method='Brent',solver=solver,rad=rad,slqd=slqd),
      'NM' = optim(tau, mll, type=type, beta=beta,u=u,s_d=s_d,sigma_s=corr,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,eps=eps,method='Nelder-Mead',solver=solver,rad=rad,slqd=slqd),
      stop("The argument opt should be bobyqa, Brent or NM.")
    )
  }else{
    new_t = switch(
      opt,
      'bobyqa' = bobyqa(tau, mll_mm, type=type, beta=beta,u=u,s_d=s_d,corr=corr,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,detap=detap,eps=eps,lower=rep(min_tau,ncm),upper=rep(max_tau,ncm),solver=solver,rad=rad,slqd=slqd),
      'L-BFGS-B' = optim(tau, mll_mm, type=type, beta=beta,u=u,s_d=s_d,corr=corr,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,detap=detap,eps=eps,lower=rep(min_tau,ncm),upper=rep(max_tau,ncm),method='L-BFGS-B',solver=solver,rad=rad,slqd=slqd),
      stop("The argument opt should be bobyqa or L-BFGS-B.")
    )
  }
  marg_ll = new_t$value
  if(opt=='bobyqa')
  {iter <- new_t$iter}else{
    iter <- new_t$counts
  }
  
  tau_e <- new_t$par
  check_tau(tau_e,min_tau,max_tau,ncm)
  
  if(ncm==1)
  {
    if(type=='dense')
    {
      re <- irls_ex(beta, u, tau_e, s_d, corr, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=detap,solver=solver)
    }else{
      re <- irls_fast_ap(beta, u, tau_e, s_d, corr, inv, X, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,solver=solver)
    }
  }else{
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
      re <- irls_ex(beta, u, 1, s_d, corr,sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=detap,solver=solver)
    }else{
      corr <- as(corr,'dgCMatrix')
      if(inv==TRUE)
      {
        corr <- Matrix::chol2inv(Matrix::chol(corr))
        corr <- as(corr,'dgCMatrix')
      }
      s_d <- as.vector(Matrix::diag(corr))
      re <- irls_fast_ap(beta, u, 1, s_d, corr,inv,X, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,solver=solver)
    }
  }
  
  if(k>0)
  {
    res_beta = as.vector(re$beta)
    res_var = diag(as.matrix(re$v11))
    HR = exp(res_beta)
    sdb = sqrt(res_var)
    p = as.vector(pchisq(res_beta^2/res_var,1,lower.tail=FALSE))
  }else{
    res_beta=HR=sdb=p=NULL
  }
  
  res <- list(beta=res_beta,HR=HR,sd_beta=sdb,p=p,tau=tau_e,iter=iter,rank=rk_cor,nsam=n,int_ll=marg_ll)
  return(res)
}

