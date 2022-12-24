
#' Estimate HRs using PPL given a known variance component (tau)
#'
#' \code{fit_ppl} returns estimates of HRs and their p-values given a known variance component (tau).
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
#' 
#' @param tau A positive scalar or vector. A variance component(s) given by the user. If there are more than one related matrix, this must be a vector, the length of which corresponds to the number of matrices. 
#' @param X A matrix of the predictors. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns.
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for failure / 0 for censored).
#' @param corr A relatedness matrix or a List object of matrices if there are multiple relatedness matrices. They can be a matrix or a 'dgCMatrix' class in the Matrix package. The matrix (or the sum if there are multiple) must be symmetric positive definite or symmetric positive semidefinite. The order of subjects must be consistent with that in outcome.
#' @param type A string indicating the sparsity structure of the relatedness matrix. Should be 'bd' (block diagonal), 'sparse', or 'dense'. See details.
#' @param eps An optional positive value indicating the relative convergence tolerance in the optimization algorithm. Default is 1e-6.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. 
#' @param solver An optional binary value that can be either 1 (Cholesky Decomposition using RcppEigen) or 2 (PCG). Default is NULL, which lets the function select a solver. See details.
#' @param verbose An optional logical value indicating whether to print additional messages. Default is TRUE.
#' @param order An optional integer value starting from 0. Only valid when dense=FALSE. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return p: The p-value.
#' @return iter: The number of iterations until convergence.
#' @return ppl: The PPL when the convergence is reached.
#' @keywords Cox mixed-effects model
#' @export fit_ppl
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
#' ## simulate random effexts and outcomes
#' x <- mvrnorm(1, rep(0,n), tau_var*sigma)
#' myrates <- exp(x-1)
#' y <- rexp(n, rate = myrates)
#' cen <- rexp(n, rate = 0.02 )
#' ycen <- pmin(y, cen)
#' outcome <- cbind(ycen,as.numeric(y <= cen))
#' 
#' ## fit the ppl
#' re = fit_ppl(x,outcome,sigma,type='bd',tau=0.5,order=1)
#' re

fit_ppl <- function(X,outcome,corr,type,tau,eps=1e-6,order=1,solver=NULL,spd=TRUE,verbose=TRUE){

  if(eps<0)
  {eps <- 1e-6}
  
  if(!(type %in% c('bd','sparse','dense')))
  {stop("The type argument should be 'bd', 'sparse' or 'dense'.")}
  
  ncm <- 1
  if(is.list(corr))
  {
    ncm <- length(corr)
    ntau <- length(tau)
    if(ncm != ntau)
    {stop("The length of tau must be the same as the number of matrices in corr.")}
    if(ncm==1)
    {
      corr <- corr[[1]]
    }else{
      for(i in 1:ncm)
      {
        corr[[i]] <- corr[[i]]*tau[i]
      }
      corr <- Reduce('+',corr)
      tau <- 1
    }
  }
  if((spd==FALSE) & (ncm>1))
  {stop("The option spd=FALSE does not support more than one correlation matrix. If multiple correlation matrices are provided, please make sure that their sum is SPD.")}
  
  X <- as.matrix(X)
  outcome <- as.matrix(outcome)
  
  if(nrow(outcome)!=nrow(X))
  {stop("The phenotype and predictor matrices have different sample sizes.")}
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem, ,drop = FALSE]
    X <- as.matrix(X[-rem,,drop = FALSE])
    corr <- corr[-rem,-rem,drop = FALSE]
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  x_sd = which(as.vector(apply(X,2,sd))>0)
  x_ind = length(x_sd)
  if(x_ind==0)
  {stop("The predictors are all constants after the removal of subjects.")}else{
    k <- ncol(X)
    if(x_ind<k)
    {
      warning(paste0(k-x_ind," predictor(s) is/are removed because they are all constants after the removal of subjects."))
      X = X[,x_sd,drop=FALSE]
      k <- ncol(X)
    }
  }
  
  n <- nrow(outcome)
    
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  u <- rep(0,n)
  beta <- rep(0,k)
  
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
  
  nz <- nnzero(corr)
  if( nz > ((as.double(n)^2)/2) )
  {type <- 'dense'}
  inv = NULL
  
  eigen = TRUE
  sigma_i_s <- NULL
  if(type=='dense')
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
    inv <- TRUE
    s_d <- as.vector(diag(sigma_i_s))
    
    solver <- get_solver(solver,type,verbose)
    
  }else{
    
    corr <- as(corr, 'dgCMatrix')
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
      corr <- Matrix::chol2inv(Matrix::chol(corr))
      inv = TRUE
    }else{
      inv = FALSE
    }
    s_d <- as.vector(Matrix::diag(corr))
    
    solver <- get_solver(solver,type,verbose)
    
    if(inv==TRUE)
    {
      corr <- as(corr,'dgCMatrix')
    }
  }
  
  if(verbose==TRUE)
  {
    model_info(inv,ncm,eigen,type,solver)
  }
  
  if(type=='dense')
  {
    res <- irls_ex(beta, u, tau, s_d, corr, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap='exact',solver=solver)
  }else{
    res <- irls_fast_ap(beta, u, tau, s_d, corr, inv, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,order,det=FALSE,detap='exact',solver=solver)
  }
  
  res_beta = as.vector(res$beta)
  res_var = diag(as.matrix(res$v11))
  p = pchisq(res_beta^2/res_var,1,lower.tail=FALSE)
  
  re = list(beta=res_beta,HR=exp(res_beta),sd_beta=sqrt(res_var),p=as.vector(p),iter=res$iter,ppl=res$ll)
  return(re)
}

