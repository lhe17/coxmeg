// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

//
// via the exports attribute we tell Rcpp to make this function
// available from R

// [[Rcpp::export]]
Eigen::MatrixXd csqei(const Eigen::Map<Eigen::VectorXd> & w_v, const Eigen::MatrixXd & mx, const Eigen::Map<Eigen::VectorXd> & rs_rs, 
                      const Eigen::Map<Eigen::VectorXd> & rs_cs,const Eigen::MatrixXi & ind, const Eigen::Map<Eigen::VectorXd> & av) {
  // multiplication of the Lap matrix (dense part) and mx
  int nr = mx.rows();
  int nc = mx.cols();
  
  Eigen::MatrixXd w_v_x = mx.array().colwise() * w_v.array();
  Eigen::MatrixXd temp(nr,nc);
  for(int j=0; j<nr; j++)
    temp.row(j) = w_v_x.row(ind(j,0));
  temp = temp.colwise().reverse().eval();
  
  for(int j=1; j<nr; ++j)
    temp.row(j) += temp.row(j-1);
  for(int j=0; j<nr; j++)
    w_v_x.row(j) = temp.row(rs_rs(j));
  w_v_x = w_v_x.array().colwise() * av.array();
  for(int j=1; j<nr; ++j)
    w_v_x.row(j) += w_v_x.row(j-1);
  for(int j=0; j<nr; j++)
    temp.row(j) = w_v_x.row(rs_cs(j));	
  for(int j=0; j<nr; j++)
    w_v_x.row(j) = temp.row(ind(j,1));
  w_v_x = w_v_x.array().colwise() * w_v.array();
  return w_v_x;
  
}



// [[Rcpp::export]]
Eigen::MatrixXd wma_cp(const Eigen::Map<Eigen::VectorXd> & w, const Eigen::Map<Eigen::VectorXd> & cs_p, const Eigen::MatrixXi & ind,
                       const Eigen::Map<Eigen::VectorXd> & a) {
  // the Lap matrix (dense part) 
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
  
  Eigen::MatrixXd q(n,n);
  for(int j=0;j<n;j++)
  {
    for(int k=j;k<n;k++)
    {
      int min = cs_p(ind(j,1));
      int ind_2 = ind(k,1);
      if(min>cs_p(ind_2))
        min = cs_p(ind_2);
      
      double q_t = w(j)*w(k)*a2csum(min);
      q(j,k) = q_t;
      q(k,j) = q_t;
    }
  }
  
  return q;
}


//
// [[Rcpp::export]]
Eigen::MatrixXd score_test(const Eigen::Map<Eigen::VectorXd> & deriv, const Eigen::Map<Eigen::VectorXd> & bw_v,
               const Eigen::Map<Eigen::VectorXd> &  w, const Eigen::Map<Eigen::VectorXd> & rs_rs,const Eigen::Map<Eigen::VectorXd> & rs_cs,
               const Eigen::Map<Eigen::VectorXd> & cs_p, const Eigen::MatrixXi & ind,
               const Eigen::Map<Eigen::VectorXd> & a, const Eigen::Map<Eigen::VectorXd> & a2, const Eigen::VectorXd & tau, 
               const Eigen::Map<Eigen::MatrixXd> & v,const Eigen::MatrixXd & cov, const Eigen::MatrixXd & x) {

    int n = bw_v.size();
    int n_c = cov.cols();
    if(cov.rows()==0)
    {
      n_c = 0;
    }
    
    int ns = x.cols();  
    Eigen::MatrixXd tst(ns,2);
    for(int i=0; i<ns; i++)
    {
      Eigen::RowVectorXd xh = bw_v.array()*x.col(i).array() - csqei(w,x.col(i),rs_rs,rs_cs,ind,a2).col(0).array(); 
      Eigen::RowVectorXd v21(n+n_c);
      if(n_c>0)
      {
        v21.head(n_c) = xh*cov;
      }
      v21.tail(n) = xh;
      double sc_v = (xh*x.col(i) - v21*v*v21.transpose()).value();
      double deriv_x = x.col(i).transpose()*deriv;
      tst(i,0) = deriv_x/sc_v;
      tst(i,1) = deriv_x*deriv_x/sc_v;
    }
    
    return tst;
}

