// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SimplicialLDLT<SpMat> SpChol;
typedef Eigen::LDLT<Eigen::MatrixXd> Chol;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double logdeth(Eigen::SparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> dv, const Eigen::Map<Eigen::VectorXd> bw_v,
                  const Eigen::Map<Eigen::VectorXd> w, const Eigen::Map<Eigen::VectorXd> cs_p, const Eigen::MatrixXi & v4,
                  const Eigen::Map<Eigen::VectorXd> a, const Eigen::VectorXd & tau, const Eigen::VectorXi & inv, const Eigen::VectorXi & detap) {
  
  int n = w.size();
  
  Eigen::ArrayXd a2 = a.array()*a.array();
  int as = a2.size();
  Eigen::ArrayXd a2csum(as);
  double cst = 0;
  for(int j = 0; j<as;j++)
  {
    cst += a2(j);
    a2csum(j) = cst;
  }
  
  if(detap(0)==1)
  {
    Eigen::ArrayXd qd(n);
    
    for(int j=0;j<n;j++)
    {
      qd(j) = w(j)*w(j)*a2csum(cs_p(v4(j,1)));
    }
   
    if(inv(0)==1)
    {
      for (int k=0; k<A.outerSize(); ++k)
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
          if(it.row()==it.col())
          {
            it.valueRef() = dv(it.row()) + (bw_v(it.row()) - qd(it.row()))*tau(0);
            break;
          }
        }
        
      SpChol solver(A);
      Eigen::ArrayXd logdet = solver.vectorD();
      logdet = logdet.log();
      return logdet.sum();
    }else{
      for (int k=0; k<A.outerSize(); ++k)
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
          if(it.row()==it.col())
          {
            it.valueRef() = dv(it.row()) + 1/(bw_v(it.row()) - qd(it.row()))/tau(0);
            break;
          }
        }
        
        SpChol solver(A);
      Eigen::ArrayXd logdet = solver.vectorD();
      logdet = logdet.log();
      Eigen::ArrayXd dd = bw_v.array() - qd;
      dd = dd.log();
      return logdet.sum() + dd.sum();
    }
    
  }else{
    Eigen::MatrixXd q(n,n);
    
    for(int j=0;j<n;j++)
    {
      for(int k=j;k<n;k++)
      {
        int min = cs_p(v4(j,1));
        int ind_2 = v4(k,1);
        if(min>cs_p(ind_2))
          min = cs_p(ind_2);
        double q_t = w(j)*w(k)*a2csum(min);
        q(j,k) = q_t;
        q(k,j) = q_t;
      }
    }
    
    q = Eigen::MatrixXd(A)/tau(0) - q;
    q.diagonal() += bw_v; 
    Chol solver(q);
    Eigen::ArrayXd logdet = solver.vectorD();
    logdet = logdet.log();
    return logdet.sum();
  }
  
}


// Different from logdeth, this function calculates two log-determinants
// [[Rcpp::export]]
double logdethmcmdense(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & dv, const Eigen::Map<Eigen::VectorXd> & bw_v,
               const Eigen::Map<Eigen::VectorXd> & w, const Eigen::Map<Eigen::VectorXd> & cs_p, const Eigen::MatrixXi & v4,
               const Eigen::Map<Eigen::VectorXd> & a) {
  
  int n = w.size();
  
  Eigen::ArrayXd a2 = a.array()*a.array();
  int as = a2.size();
  Eigen::ArrayXd a2csum(as);
  double cst = 0;
  for(int j = 0; j<as;j++)
  {
    cst += a2(j);
    a2csum(j) = cst;
  }
  
  Chol solver;
  solver.compute(A);
  Eigen::ArrayXd logdet2 = solver.vectorD();
  logdet2 = logdet2.log();
  
  Eigen::MatrixXd q(n,n);
  
  for(int j=0;j<n;j++)
  {
    for(int k=j;k<n;k++)
    {
      int min = cs_p(v4(j,1));
      int ind_2 = v4(k,1);
      if(min>cs_p(ind_2))
        min = cs_p(ind_2);
      double q_t = w(j)*w(k)*a2csum(min);
      q(j,k) = q_t;
      q(k,j) = q_t;
    }
  }
  
  q = A - q;
  q.diagonal() += bw_v; 
  solver.compute(q);
  Eigen::ArrayXd logdet = solver.vectorD();
  logdet = logdet.log();
  return (logdet.sum() - logdet2.sum());
  
}


// [[Rcpp::export]]
double logdet_ch(const Eigen::Map<Eigen::MatrixXd> & X_m, const Eigen::Map<Eigen::MatrixXd> & rad_m, const Eigen::Map<Eigen::VectorXd> & bma_d,
                 const Eigen::Map<Eigen::VectorXd> & bpa_d, const Eigen::Map<Eigen::VectorXd> & cj_v) {
  
  int t = rad_m.cols();
  int q = cj_v.size();
  
  double app = 0;
  double ratio = bpa_d[0]/bma_d[0];
  
  Eigen::MatrixXd u = cj_v[0]*rad_m;
  Eigen::MatrixXd w0 = rad_m;
  Eigen::MatrixXd w1 = 2/bma_d[0]*(X_m*rad_m) - ratio*rad_m;
  Eigen::MatrixXd w2 = w1;
  u = u+cj_v[1]*w1;
  for(int j=2; j<q; j++)
  {
    w2 = 4/bma_d[0]*(X_m*w1) - w0 - 2*ratio*w1;
    u = u + cj_v[j]*w2;
    w0 = w1;
    w1 = w2;
  }
  u = rad_m.array()*u.array();
  app = u.sum()/t;
  
  return app;
}

// [[Rcpp::export]]
double logdet_lanczos(const Eigen::Map<Eigen::MatrixXd> & X_m, const Eigen::Map<Eigen::MatrixXd> & rad_m, const Eigen::VectorXi & m_d) {
  
  // SLQ for a dense matrix
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::ArrayXXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m-1,t);
  Eigen::ArrayXXd temp(n,t);
  Eigen::MatrixXd T;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = X_m*v1;
  temp = w*v1.array();
  alpha.row(0) = temp.colwise().sum();
  w = w - v1.array().rowwise()*alpha.row(0);
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    temp = w*w;
    beta.row(j-1) = sqrt(temp.colwise().sum());
    
    for(int k=0;k<t;k++)
    {
      if(beta(j-1,k)<1e-10)
      {
        if(iter(k)==m)
          iter(k) = j;
      }
    }
    
    v2 = w.rowwise()*beta.row(j-1).inverse();
    
    w = X_m*v2;
    temp = w*v2.array();
    alpha.row(j) = temp.colwise().sum();
    w = w - v2.array().rowwise()*alpha.row(j) - v1.array().rowwise()*beta.row(j-1);
    v1 = v2;
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    T = Eigen::MatrixXd::Zero(it,it);
    T.diagonal() = alpha.col(i).head(it);
    if(it>1)
    {
      T.diagonal(1) = beta.col(i).head(it-1);
      T.diagonal(-1) = beta.col(i).head(it-1);
    }
    es.compute(T);
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    eivl = eivt*eivl;
    logdet = logdet + eivl.sum();
  }
  
  return logdet*n/t;
}

// [[Rcpp::export]]
Eigen::ArrayXXd wma_mv(const Eigen::MatrixXd & X_m, const Eigen::Map<Eigen::VectorXd> & w_v, const Eigen::Map<Eigen::VectorXd> & rs_rs, 
                       const Eigen::Map<Eigen::VectorXd> & rs_cs,const Eigen::Map<Eigen::MatrixXd> & ind, const Eigen::Map<Eigen::VectorXd> & av2,
                       const int n, const int t) 
{
  // almost same as csqei, just return an Array
  Eigen::MatrixXd temp(n,t);
  Eigen::ArrayXXd w_v_x(n,t);
  
  w_v_x = X_m.array().colwise() * w_v.array();
  for(int j=0; j<n; j++)
    temp.row(j) = w_v_x.row(ind(j,0));
  temp = temp.colwise().reverse().eval();
  for(int j=1; j<n; ++j)
    temp.row(j) += temp.row(j-1);
  for(int j=0; j<n; j++)
    w_v_x.row(j) = temp.row(rs_rs(j));
  w_v_x = w_v_x.colwise() * av2.array();
  for(int j=1; j<n; ++j)
    w_v_x.row(j) += w_v_x.row(j-1);
  for(int j=0; j<n; j++)
    temp.row(j) = w_v_x.row(rs_cs(j));	
  for(int j=0; j<n; j++)
    w_v_x.row(j) = temp.row(ind(j,1));
  w_v_x = w_v_x.colwise() * w_v.array();
  return w_v_x;
}

// [[Rcpp::export]]
double logdet_gkb(const Eigen::Map<Eigen::MatrixXd> & X_m, const Eigen::Map<Eigen::VectorXd> & bw_v, const double & tau, 
                  const Eigen::Map<Eigen::VectorXd> & w_v, const Eigen::Map<Eigen::VectorXd> & rs_rs, 
                  const Eigen::Map<Eigen::VectorXd> & rs_cs,const Eigen::Map<Eigen::MatrixXd> & ind, const Eigen::Map<Eigen::VectorXd> & av2,
                  const Eigen::Map<Eigen::MatrixXd> & rad_m, const Eigen::VectorXi & m_d) {
  
  // GKB for a dense matrix
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::ArrayXXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m,t);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = v1.array().colwise()*bw_v.array() - wma_mv(v1,w_v,rs_rs,rs_cs,ind,av2,n,t);
  w = X_m*w.matrix();
  if(tau!=1)
  {
    w = w + v1.array()/tau;
  }else{
    w = w + v1.array();
  }
  alpha.row(0) = (w.cwiseProduct(w)).colwise().sum();
  alpha.row(0) = alpha.row(0).array().sqrt();
  v2 = w.rowwise()*alpha.row(0).inverse();
  
  w = X_m*v2;
  w = w.colwise()*bw_v.array() - wma_mv(w,w_v,rs_rs,rs_cs,ind,av2,n,t);
  if(tau!=1)
  {
    w = w + v2.array()/tau;
  }else{
    w = w + v2.array();
  }
  w = w - v1.array().rowwise()*alpha.row(0);
  
  beta.row(0) = (w.cwiseProduct(w)).colwise().sum();
  beta.row(0) = beta.row(0).array().sqrt();
  v1 = w.rowwise()*beta.row(0).inverse();
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    w = v1.array().colwise()*bw_v.array() - wma_mv(v1,w_v,rs_rs,rs_cs,ind,av2,n,t);
    w = X_m*w.matrix();
    if(tau!=1)
    {
      w = w + v1.array()/tau;
    }else{
      w = w + v1.array();
    }
    w = w - v2.array().rowwise()*beta.row(j-1);
    
    alpha.row(j) = (w.cwiseProduct(w)).colwise().sum();
    alpha.row(j) = alpha.row(j).array().sqrt();
    v2 = w.rowwise()*alpha.row(j).inverse();
    
    w = X_m*v2;
    w = w.colwise()*bw_v.array() - wma_mv(w,w_v,rs_rs,rs_cs,ind,av2,n,t);
    if(tau!=1)
    {
      w = w + v2.array()/tau;
    }else{
      w = w + v2.array();
    }
    w = w - v1.array().rowwise()*alpha.row(j);
    
    beta.row(j) = (w.cwiseProduct(w)).colwise().sum();
    beta.row(j) = beta.row(j).array().sqrt();
    v1 = w.rowwise()*beta.row(j).inverse();
    
    for(int k=0;k<t;k++)
    {
      if((beta(j,k)<1e-10) || (alpha(j,k)<1e-10))
      {
        if(iter(k)==m)
          iter(k) = j+1;
      }
    }
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    Eigen::VectorXd maindi = alpha.col(i).head(it).pow(2);
    maindi.tail(it-1) = maindi.tail(it-1).array() + beta.col(i).head(it-1).pow(2);
    Eigen::VectorXd subdi = alpha.col(i).head(it-1)*beta.col(i).head(it-1);
    es.computeFromTridiagonal(maindi,subdi);
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    logdet = logdet + (eivt*eivl).sum();
  }
  
  return logdet*n/t;
}

// [[Rcpp::export]]
double logdet_lanczos_sp(const Eigen::MappedSparseMatrix<double> & X_m, const Eigen::Map<Eigen::MatrixXd> & rad_m, const Eigen::VectorXi & m_d) {
  
  // SLQ for a sparse matrix
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::ArrayXXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m-1,t);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = X_m*v1;
  alpha.row(0) = (w*v1.array()).colwise().sum();
  w = w - v1.array().rowwise()*alpha.row(0);
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    beta.row(j-1) = (w*w).colwise().sum();
    beta.row(j-1) = beta.row(j-1).sqrt();
    
    for(int k=0;k<t;k++)
    {
      if(beta(j-1,k)<1e-10)
      {
        if(iter(k)==m)
          iter(k) = j;
      }
    }
    
    v2 = w.rowwise()*beta.row(j-1).inverse();
    
    w = X_m*v2;
    alpha.row(j) = (w*v2.array()).colwise().sum();
    w = w - v2.array().rowwise()*alpha.row(j) - v1.array().rowwise()*beta.row(j-1);
    v1 = v2;
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    es.computeFromTridiagonal(alpha.col(i).head(it),beta.col(i).head(it-1));
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    eivl = eivt*eivl;
    logdet = logdet + eivl.sum();
  }
  
  return logdet*n/t;
}


// [[Rcpp::export]]
double logdet_gkb_sp(const Eigen::MappedSparseMatrix<double> & X_m, const Eigen::Map<Eigen::MatrixXd> & rad_m, const Eigen::VectorXi & m_d) {
  
  // GKB for a sparse matrix
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::MatrixXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m,t);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = X_m*v1;
  alpha.row(0) = (w.cwiseProduct(w)).colwise().sum();
  alpha.row(0) = alpha.row(0).array().sqrt();
  v2 = w.array().rowwise()*alpha.row(0).inverse();
  w = (X_m.transpose()*v2).array() - v1.array().rowwise()*alpha.row(0);

  beta.row(0) = (w.cwiseProduct(w)).colwise().sum();
  beta.row(0) = beta.row(0).array().sqrt();
  v1 = w.array().rowwise()*beta.row(0).inverse();
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    w = (X_m*v1).array() - v2.array().rowwise()*beta.row(j-1);
    
    alpha.row(j) = (w.cwiseProduct(w)).colwise().sum();
    alpha.row(j) = alpha.row(j).array().sqrt();
    v2 = w.array().rowwise()*alpha.row(j).inverse();
    
    w = (X_m.transpose()*v2).array() - v1.array().rowwise()*alpha.row(j);

    beta.row(j) = (w.cwiseProduct(w)).colwise().sum();
    beta.row(j) = beta.row(j).array().sqrt();
    v1 = w.array().rowwise()*beta.row(j).inverse();
    
    for(int k=0;k<t;k++)
    {
      if((beta(j,k)<1e-10) || (alpha(j,k)<1e-10))
      {
        if(iter(k)==m)
          iter(k) = j+1;
      }
    }
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    Eigen::VectorXd maindi = alpha.col(i).head(it).pow(2);
    maindi.tail(it-1) = maindi.tail(it-1).array() + beta.col(i).head(it-1).pow(2);
    Eigen::VectorXd subdi = alpha.col(i).head(it-1)*beta.col(i).head(it-1);
    es.computeFromTridiagonal(maindi,subdi);
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    logdet = logdet + (eivt*eivl).sum();
  }
  
  return logdet*n/t;
}

// [[Rcpp::export]]
double logdet_gkb_dense(const Eigen::Map<Eigen::MatrixXd> & X_m, const Eigen::Map<Eigen::MatrixXd> & rad_m, const Eigen::VectorXi & m_d) {
  
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::MatrixXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m,t);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = X_m*v1;
  alpha.row(0) = (w.cwiseProduct(w)).colwise().sum();
  alpha.row(0) = alpha.row(0).array().sqrt();
  v2 = w.array().rowwise()*alpha.row(0).inverse();
  w = (X_m.transpose()*v2).array() - v1.array().rowwise()*alpha.row(0);
  
  beta.row(0) = (w.cwiseProduct(w)).colwise().sum();
  beta.row(0) = beta.row(0).array().sqrt();
  v1 = w.array().rowwise()*beta.row(0).inverse();
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    w = (X_m*v1).array() - v2.array().rowwise()*beta.row(j-1);
    
    alpha.row(j) = (w.cwiseProduct(w)).colwise().sum();
    alpha.row(j) = alpha.row(j).array().sqrt();
    v2 = w.array().rowwise()*alpha.row(j).inverse();
    
    w = (X_m.transpose()*v2).array() - v1.array().rowwise()*alpha.row(j);
    
    beta.row(j) = (w.cwiseProduct(w)).colwise().sum();
    beta.row(j) = beta.row(j).array().sqrt();
    v1 = w.array().rowwise()*beta.row(j).inverse();
    
    for(int k=0;k<t;k++)
    {
      if((beta(j,k)<1e-10) || (alpha(j,k)<1e-10))
      {
        if(iter(k)==m)
          iter(k) = j+1;
      }
    }
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    Eigen::VectorXd maindi = alpha.col(i).head(it).pow(2);
    maindi.tail(it-1) = maindi.tail(it-1).array() + beta.col(i).head(it-1).pow(2);
    Eigen::VectorXd subdi = alpha.col(i).head(it-1)*beta.col(i).head(it-1);
    es.computeFromTridiagonal(maindi,subdi);
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    logdet = logdet + (eivt*eivl).sum();
  }
  
  return logdet*n/t;
}

