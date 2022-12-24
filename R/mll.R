mll <- function(tau, type, beta, u,s_d,sigma_s, sigma_i_s,X, d_v, ind, rs_rs, rs_cs,rs_cs_p,det,detap,inv,rk_cor,eps=1e-6,order=0,solver=1,rad=NULL,slqd=8){
  
  if(type=='dense')
  {
    re <- irls_ex(beta, u, tau, s_d, sigma_s, sigma_i_s, X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,det=det,detap=detap,solver=solver,rad=rad,slqd=slqd)
  }else{
      re <- irls_fast_ap(beta, u, tau, s_d, sigma_s, inv, X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,order,det=det,detap=detap,solver=solver,rad=rad,slqd=slqd)
  }

  return(as.numeric(rk_cor*log(tau) + re$logdet[1] + (-2)*re$ll))
}


mll_mm <- function(tau, type, beta, u,s_d,corr, X, d_v, ind, rs_rs, rs_cs,rs_cs_p,detap,rk_cor,eps=1e-6,order=0,solver=1,rad=NULL,slqd=8){
  
  ncm <- length(corr)
  for(i in 1:ncm)
  {
    corr[[i]] <- corr[[i]]*tau[i]
  }
  corr <- as.matrix(Reduce('+',corr))
  
  if(type=='dense')
  {
    sigma_i_s = chol(corr)
    sigma_i_s = as.matrix(chol2inv(sigma_i_s))
    s_d <- as.vector(diag(sigma_i_s))
    re <- irls_ex(beta, u, 1, s_d, corr, sigma_i_s, X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,det=TRUE,detap=detap,solver=solver,rad=rad,slqd=slqd)
  }else{
    corr <- as(corr,'dgCMatrix')
    s_d <- as.vector(Matrix::diag(corr))
    re <- irls_fast_ap(beta, u, 1, s_d, corr, FALSE, X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,order,det=TRUE,detap=detap,solver=solver,rad=rad,slqd=slqd)
  }
  
  return(as.numeric(re$logdet[1] + (-2)*re$ll))
}


check_input <- function(outcome,corr,type,detap,ncm,spd)
{
  if(!(type %in% c('bd','sparse','dense')))
  {stop("The type argument should be 'bd', 'sparse' or 'dense'.")}
  
  if(!is.null(detap))
  {
    if(!(detap %in% c('exact','slq','diagonal','gkb')))
    {stop("The detap argument should be 'exact', 'diagonal', 'gkb' or 'slq'.")}
  }
  
  if((spd==FALSE) & (ncm>1))
  {stop("The option spd=FALSE does not support more than one correlation matrix. If multiple correlation matrices are provided, please make sure that their sum is SPD.")}
  
  nro = nrow(outcome)
  if(ncm==1)
  {
    if((nro!=nrow(corr)) || (nro!=ncol(corr)))
    {stop("The phenotype and the relatedness matrix have different sample sizes.")}
  }else{
    nrowc = unlist(lapply(corr,function(x)nrow(x)))
    ncolc = unlist(lapply(corr,function(x)ncol(x)))
    if(length(unique(c(nro,nrowc,ncolc)))>1)
    {
      stop("The phenotype and all relatedness matrices must have the same sample size.")
    }
  }
  
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  
}

model_info <- function(inv,ncm,eigen,type,solver,detap=NULL)
{
  if(type=='dense')
  {
    message('The correlation matrix is treated as dense.')
  }else{
    message('The correlation matrix is treated as sparse/block diagonal.')
  }
  
  if((inv==TRUE) & (ncm==1))
  {message('The relatedness matrix is inverted.')}
  
  if(is.null(detap)==FALSE)
  {
    if(detap=='mcm')
    {detap <- 'exact'}
    message("The method for computing the determinant is '", detap, "'.")
  }
  
  if(type=='dense')
  {
    switch(
      solver,
      '1' = message('Solver: solve (base).'),
      '2' = message('Solver: PCG (RcppEigen:dense).')
    )
  }else{
    switch(
      solver,
      '1' = message('Solver: Cholesky decomposition (RcppEigen=',eigen,').'),
      '2' = message('Solver: PCG (RcppEigen:sparse).')
    )
  }
}

get_solver <- function(solver,type,verbose)
{
  if(type=='dense')
  {
    if(is.null(solver))
    {solver = 2}else{
      if(!(solver%in%c(1,2)))
      {
        solver <- 2
        if(verbose==TRUE)
        {message("Invalid value for solver. solver=2 is used.")}
      }
    }
  }else{
    if(is.null(solver))
    {
      if(type=='bd')
      {
        solver = 1
      }else{solver = 2}
    }else{
      if(!(solver%in%c(1,2)))
      {
        solver <- 2
        if(verbose==TRUE)
        {message("Invalid value for solver. solver=2 is used.")}
      }
    }
  }
  solver
}

set_detap_dense <- function(detap,n,spsd,ncm)
{
  if(is.null(detap))
  {
    if(ncm==1)
    {
      if(n<3000)
      {detap = 'exact'}else{
        #if(spsd==FALSE)
        #{detap = 'slq'}else{detap = 'gkb'}
        detap = 'gkb'
      }
    }else{
      if(n<3000)
      {detap = 'mcm'}else{detap = 'gkb'}
    }
  }else{
    if(ncm>1)
    {
      if(detap=='exact')
      {detap <- 'mcm'}else{detap <- 'gkb'}
    }
  }
  detap
}

set_detap_sparse <- function(detap,type,verbose,ncm)
{
  if(is.null(detap))
  {
    detap <- ifelse(type=='bd','diagonal','slq')
  }else{
    if((detap=='exact') & (ncm>1 | type=='sparse'))
    {
      detap = 'slq'
      if(verbose==TRUE)
      {message("detap=exact is not supported under this setting. The detap argument is changed to 'slq'.")}
    }
  }
  detap
}

check_tau <- function(tau_e,min_tau,max_tau,ncm)
{
  if(sum(tau_e==min_tau)==ncm)
  {warning(paste0("The estimated variance component equals the lower bound (", min_tau, "), probably suggesting no random effects."))}
  if(sum(tau_e==max_tau)>0)
  {warning(paste0("At least one of the estimated variance component equals the upper bound (", max_tau, "). It might be better to use a larger max_tau."))}
  
}
