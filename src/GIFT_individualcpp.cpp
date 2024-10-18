#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


void lmm_pxem_ptr2_individual(const arma::vec& X, const arma::mat& Zx, arma::mat Omega,  const int& maxIter,
                     arma::vec& SigGdiag, arma::vec& d, const arma::vec& constrp,int k, double& loglik_max,
                     int& iteration, arma::mat& Sigb, arma::vec& mub){
  
  int n = X.n_elem/k, p = Zx.n_cols ;
  
  if (p != (int)mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }
  
  if (p != (int)Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }
    
  double lambda, lambda2;  // parameter expansion
  
  // initialize
  vec contrn=zeros<vec>(k+1);
  int su=0;
  for(int v=0; v < k+1; v++){	  
    contrn(v) = su;
    su += n;
  }
  
  for(int v=0; v < k; v++){
    SigGdiag(v) = 1.0;
  }

  int r=0;
  vec rr(k);
  for(int v=0; v < k; v++){
	double dd = SigGdiag(v)/constrp(v);
	rr(v) = dd;
    for(int h=0; h < constrp(v); h++){
	  d(r) += dd;
	  r++;
    }
  }

  vec constrp1(k+1);
  constrp1(0)= 0;
  for (int v = 0; v < k; v++){
	constrp1(v+1)=constrp(v);
  }
  
  int sum0=0;
  vec constrp2(k+1);
  constrp2(0)= 0;
  for (int v = 0; v < k; v++){
	sum0 += constrp(v);
	constrp2(v+1) = sum0;
  }

  vec loglik(maxIter);
  loglik(0) = -datum::inf;
  
  mat Sigb_inv(p,p),Op = ones<mat>(p,1), K1(p,p), Omega_inv(k,k), OR(k,k),invOR(k,k);

  double  E, M;

  OR = chol(Omega);
  invOR = inv(OR);
  Omega_inv = invOR*invOR.t(); 

  mat XK23=zeros<mat>(1,p+1);
  for(int v=0; v < k; v++){	 
	mat XK21=zeros<mat>(constrp(v),1);
    if(v!=0){
	  XK21 = zeros<mat>(constrp(v),constrp2(v)+1);
    }
    for(int u=v; u < k; u++){	 
        mat XK22 = Zx.submat(0,constrp2(v),n-1,constrp2(v+1)-1).t()*Zx.submat(0,constrp2(u),n-1,constrp2(u+1)-1) * as_scalar(Omega_inv(v,u));
        XK21 = join_rows(XK21,XK22); 
    }
    XK23 =join_cols(XK23,XK21);
  }
  XK23.shed_col(0);
  XK23.shed_row(0);
  mat XK24 = symmatu(XK23);

  mat KX23=zeros<mat>(1,k+1);
  for(int v=0; v < k; v++){	 
  mat KX21=zeros<mat>(constrp(v),1);
  for(int u=0; u < k; u++){	 
  mat KX22 = Zx.submat(0,constrp2(v),n-1,constrp2(v+1)-1).t()*X.subvec(contrn(u),contrn(u+1)-1);
  KX21 = join_rows(KX21,KX22); 
  }
  KX23 =join_cols(KX23,KX21);
  }
  KX23.shed_col(0);
  KX23.shed_row(0);

  mat RR(p,p),invRR(p,p);

  iteration = maxIter-1;
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
	mat B1=zeros(1,1);
	vec O2=1.0/SigGdiag;
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=as_scalar(O2(v))*O1;
	  B1=join_cols(B1,O);
    }
    B1.shed_row(0);
 
    K1 = (B1*B1.t()) % XK24;

    Sigb_inv = K1;
	Sigb_inv.diag()+= 1.0/d;

    RR = chol(Sigb_inv);
    invRR = inv(RR);
    Sigb = invRR*invRR.t();

	mat AA = Omega_inv % (O2*O2.t());
	mat tmp = KX23 * AA;
    mat in = tmp.submat(0,0,constrp(0)-1,0);
	if(k>1){
    for(int u=1; u < k; u++){	 
        mat XK2 = tmp.submat(constrp2(u),u,constrp2(u+1)-1,u);
        in = join_cols(in,XK2); 
    }
    }
    mub = Sigb*in;

    E = (n-1) * accu(Omega % Omega_inv % (O2*O2.t())) - 2*sum(mub % in) + sum(mub%(K1*mub));
    M = accu(mub % (1.0/d) % mub);
	
	loglik(iter) = - n*sum(log(SigGdiag % SigGdiag))*0.5 - sum(log(d))/2 - E/2- M/2 - sum(log(RR.diag())) - n*k/2*log(2*datum::pi);

    //if ( loglik(iter) - loglik(iter - 1) < 0 ){
    //  perror("The likelihood failed to increase!");
    //}
    
    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }
    
    // M-step
    lambda = sum(mub % in) / (sum(mub % (K1*mub)) + trace(K1*Sigb));
    lambda2 = pow(lambda , 2);

	for(int g = 0; g < k; g++){
	  vec G = X.subvec(contrn(g),contrn(g+1)-1)-lambda*Zx.submat(0,constrp2(g),n-1,constrp2(g+1)-1)*mub.subvec(constrp2(g),constrp2(g+1)-1);
	  double p1=sum(G%G)* as_scalar(Omega_inv(g,g)) +lambda2*trace((XK24.submat(constrp2(g),constrp2(g),constrp2(g+1)-1,constrp2(g+1)-1))*(Sigb.submat(constrp2(g),constrp2(g),constrp2(g+1)-1,constrp2(g+1)-1)));
	  double p2=0;
	  for(int u = 0; u < k; u++){
	  vec G1 =  X.subvec(contrn(u),contrn(u+1)-1)-lambda*Zx.submat(0,constrp2(u),n-1,constrp2(u+1)-1)*mub.subvec(constrp2(u),constrp2(u+1)-1);
	  p2 += ((sum(G%G1)* as_scalar(Omega_inv(g,u)) +lambda2*trace((XK24.submat(constrp2(u),constrp2(g),constrp2(u+1)-1,constrp2(g+1)-1))*(Sigb.submat(constrp2(g),constrp2(u),constrp2(g+1)-1,constrp2(u+1)-1)))))/SigGdiag(u);
	  }
	  SigGdiag(g)= 2*p1/(-(p2-p1/SigGdiag(g))+sqrt(pow((p2-p1/SigGdiag(g)),2)+4*n*p1));
	}
	
    r=0;
	for(int v=0; v < k; v++){
	  double dd=sum(mub.subvec(constrp2(v),constrp2(v+1)-1) % mub.subvec(constrp2(v),constrp2(v+1)-1))/constrp1(v+1)+ sum(Sigb.submat(constrp2(v),constrp2(v),constrp2(v+1)-1,constrp2(v+1)-1).diag())/constrp1(v+1);
	  rr(v) = dd;
      for(int h=0; h < constrp(v); h++){
	    d(r) = dd;
	    r++;
      }
    }
	
    // Reduction step
	d = lambda2 * d;

    // lambda = 1;
    // lambda2 = pow(lambda , 2);
  }
   
  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);
  
  loglik_max = loglik(iteration);
}


void loglike_twas_individual(const arma::mat& R, double& E, double& res_ysum, double& M, const arma::vec& eVal, const arma::vec& d,
                    double& sigma2y,const int& n1, const int& n2,  double& loglik){
  double Eab;//, loglik;
  Eab =E/2 +res_ysum/2/sigma2y + M/2;
  
  loglik = - log(sigma2y)*n2*0.5- accu(log(d))*0.5 - accu(log(eVal))*n1*0.5 - sum(log(R.diag())) - Eab;

}


// [[Rcpp::export]]
List GIFT_individualcpp(arma::vec X, arma::vec y, arma::mat Zx, arma::mat Zy, arma::mat Omega, arma::vec constralpha, arma::vec constrp, int maxIter, double epsStopLogLik){// *
  
  int pxem_indicator = 1;
   	
  // *convert data format
  int k = constrp.n_elem;
  int n2 = y.n_elem, p1 = Zx.n_cols, p2 = Zy.n_cols, n1 = X.n_elem/k;
  if (p1 != p2){
    perror("The dimensions of Zx and Zy are not matched");
  }
  if (k != Omega.n_cols){
    perror("The dimensions of X and Omega are not matched");
  }
  int p = p1;

    // initialization using lmm_pxem
  double sigma2y,loglik0;
  int iter0;
  mat Sigb = zeros<mat>(p,p);
  vec mub  = zeros<vec>(p), d = zeros<vec>(p), SigGdiag = zeros<vec>(k);

  lmm_pxem_ptr2_individual(X, Zx, Omega, maxIter , SigGdiag , d , constrp, k, loglik0 , iter0 , Sigb , mub);

  //declaration of variables used within loop
  mat Omega_inv(k,k), OR(k,k), invOR(k,k),Sigb_inv(p,p), R(p,p), invR(p,p), K1(p,p), K2(p,p), In1 = eye(n1,n1);
  vec mutmp(p),  mu(p),  res_y(n2), G(n1);
	
  //initialization of parameters
  double lambda = 1, lambda2 = lambda*lambda;
  mat alpha = zeros<mat>(k,1), alpha2(k,k);
  sigma2y = var(y);

  vec loglik(maxIter), loglik2(maxIter);

  int sum0=0;
  vec constrp1(k+1);
  constrp1(0)= 0;
  for (int v = 0; v < k; v++){
	sum0 += constrp(v);
	constrp1(v+1) = sum0;
  }
	  
  vec contrn=zeros<vec>(k+1);
  int su=0;
  for(int v=0; v < k+1; v++){	  
    contrn(v) = su;
    su += n1;
  }

  OR = chol(Omega);
  invOR = inv(OR);
  Omega_inv = invOR*invOR.t();

  mat K23=zeros<mat>(1,p+1);
  for(int v=0; v < k; v++){	 
	mat K21=zeros<mat>(constrp(v),1);
    if(v!=0){
	  K21 = zeros<mat>(constrp(v),constrp1(v)+1);
    }
    for(int u=v; u < k; u++){	 
      mat K22 = Zy.submat(0,constrp1(v),n2-1,constrp1(v+1)-1).t()*Zy.submat(0,constrp1(u),n2-1,constrp1(u+1)-1);
      K21 = join_rows(K21,K22); 
    }
    K23 =join_cols(K23,K21);
  }
  K23.shed_col(0);
  K23.shed_row(0);
  mat K24 = symmatu(K23);

  mat Ky23=zeros<mat>(1,1);
  for(int v=0; v < k; v++){	 
    mat K22 = Zy.submat(0,constrp1(v),n2-1,constrp1(v+1)-1).t()*y;
    Ky23 = join_cols(Ky23,K22); 
  }
  Ky23.shed_row(0);
	
  mat XK23=zeros<mat>(1,p+1);
  for(int v=0; v < k; v++){	 
    mat XK21=zeros<mat>(constrp(v),1);
    if(v!=0){
	  XK21 = zeros<mat>(constrp(v),constrp1(v)+1);
    }
    for(int u=v; u < k; u++){	 
        mat XK22 = Zx.submat(0,constrp1(v),n1-1,constrp1(v+1)-1).t()*Zx.submat(0,constrp1(u),n1-1,constrp1(u+1)-1) * as_scalar(Omega_inv(v,u));
        XK21 = join_rows(XK21,XK22); 
    }
    XK23 =join_cols(XK23,XK21);
  }
  XK23.shed_col(0);
  XK23.shed_row(0);
  mat XK24 = symmatu(XK23);
  
  mat KX23=zeros<mat>(1,k+1);
  for(int v=0; v < k; v++){	 
  mat KX21=zeros<mat>(constrp(v),1);
  for(int u=0; u < k; u++){	 
  mat KX22 = Zx.submat(0,constrp1(v),n1-1,constrp1(v+1)-1).t()*X.subvec(contrn(u),contrn(u+1)-1);
  KX21 = join_rows(KX21,KX22); 
  }
  KX23 =join_cols(KX23,KX21);
  }
  KX23.shed_col(0);
  KX23.shed_row(0);


  // initialization of likelihood
  loglik(0) = NAN;
  int Iteration = 1;
  for (int iter = 2; iter <= maxIter; iter ++ ) {
    // E-step
    vec eVal = SigGdiag % SigGdiag;
	mat B1=zeros(1,1);
	vec O2=1.0/SigGdiag;
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=as_scalar(O2(v))*O1;
	  B1=join_cols(B1,O);
    }
    B1.shed_row(0);

    K1 = (B1*B1.t()) % XK24;

	mat A1=zeros(1,1);
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=alpha(v)*O1;
	  A1=join_cols(A1,O);
    }
    A1.shed_row(0);

	K2 = ((A1*A1.t()) % K24)/sigma2y ;

	Sigb_inv = K1 + K2;
    Sigb_inv.diag() += 1.0/d;

    R = chol(Sigb_inv);
    invR = inv(R);
    Sigb = invR*invR.t();

	mat AA = Omega_inv % (O2*O2.t());
	mat tm = KX23 * AA;
    mat in = tm.submat(0,0,constrp(0)-1,0);
    if(k>1){
    for(int u=1; u < k; u++){	 
        mat XK22 = tm.submat(constrp1(u),u,constrp1(u+1)-1,u);
        in = join_cols(in,XK22); 
    }
    }
    mutmp = in + (A1 % Ky23)/sigma2y;
    mu = Sigb*mutmp;

    double E,M,res_ysum;
      
	res_ysum = n2 -1 - 2*sum(mu % (A1 % Ky23)) + as_scalar(mu.t()*K2*mu)*sigma2y;
	
    E = (n1-1) * accu(Omega % Omega_inv % (O2*O2.t())) - 2*sum(mu % in) + sum(mu%(K1*mu));
	M = accu(mu % (1.0/d)% mu);
      
    loglike_twas_individual(R, E, res_ysum, M, eVal, d, sigma2y, n1 , n2, loglik(iter - 1));
     
    //if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
    //  perror("The likelihood failed to increase!");
    //}
    // M-step
	  
    double tr1 = trace(K1*Sigb);
	   
    if (pxem_indicator == 1){
      lambda = sum(mu % in)/(sum(mu % (K1*mu)) + tr1);
      lambda2 = lambda*lambda;
    }
    else {
      lambda = 1; lambda2 = 1;
    }
	  
    uvec q =  find(constralpha == 0);
	double trde=0;
	double tmp=0;
	mat A = mu*mu.t()+Sigb;
	  
    for(int i = 0;i < k - accu(constralpha);i++){
      int a = q(i);
	  double tmp1=0;
	  for(int j = 0;j < k; j++){
		if(j != a){
          tmp1 += as_scalar(alpha(j)) * trace((K24.submat(constrp1(a),constrp1(j),constrp1(a+1)-1,constrp1(j+1)-1)) * (A.submat(constrp1(j),constrp1(a),constrp1(j+1)-1,constrp1(a+1)-1)));
		}
	  }
    trde = trace((K24.submat(constrp1(a),constrp1(a),constrp1(a+1)-1,constrp1(a+1)-1))*(A.submat(constrp1(a),constrp1(a),constrp1(a+1)-1,constrp1(a+1)-1)));
    tmp = sum(mu.subvec(constrp1(a),constrp1(a+1)-1) % (Ky23.submat(constrp1(a),0,constrp1(a+1)-1,0)));
    alpha(a) = (tmp-tmp1)/trde;
    }
    alpha2 = alpha*alpha.t();

	for(int g = 0; g < k; g++){
	  vec G = X.subvec(contrn(g),contrn(g+1)-1)-lambda*Zx.submat(0,constrp1(g),n1-1,constrp1(g+1)-1)*mu.subvec(constrp1(g),constrp1(g+1)-1);
	  double p1=(sum(G%G)* as_scalar(Omega_inv(g,g)) +lambda2*trace((XK24.submat(constrp1(g),constrp1(g),constrp1(g+1)-1,constrp1(g+1)-1))*(Sigb.submat(constrp1(g),constrp1(g),constrp1(g+1)-1,constrp1(g+1)-1))));
	  double p2=0;
	  for(int u = 0; u < k; u++){
	  vec G1 =  X.subvec(contrn(u),contrn(u+1)-1)-lambda*Zx.submat(0,constrp1(u),n1-1,constrp1(u+1)-1)*mu.subvec(constrp1(u),constrp1(u+1)-1);
	  p2 += ((sum(G%G1)* as_scalar(Omega_inv(g,u)) +lambda2*trace((XK24.submat(constrp1(u),constrp1(g),constrp1(u+1)-1,constrp1(g+1)-1))*(Sigb.submat(constrp1(g),constrp1(u),constrp1(g+1)-1,constrp1(u+1)-1)))))/SigGdiag(u);
	  }
	  SigGdiag(g)= 2*p1/(-(p2-p1/SigGdiag(g))+sqrt(pow((p2-p1/SigGdiag(g)),2)+4*n1*p1));
	}
	
	A1=zeros(1,1);
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=alpha(v)*O1;
	  A1=join_cols(A1,O);
	}
	A1.shed_row(0);
    double tr2 = trace(Sigb*((A1*A1.t()) % K24));

    res_ysum = n2 - 1 - 2*sum(mu % (A1 % Ky23)) + as_scalar(mu.t()*((A1*A1.t()) % K24)*mu);
    sigma2y = (res_ysum + tr2)/n2;

	vec rr(k);
    int r=0;
	double dd;

	for(int v=0; v < k; v++){
	  dd=sum(mu.subvec(constrp1(v),constrp1(v+1)-1) % mu.subvec(constrp1(v),constrp1(v+1)-1))/constrp(v)+ sum(Sigb.submat(constrp1(v),constrp1(v),constrp1(v+1)-1,constrp1(v+1)-1).diag())/constrp(v);
	  rr(v) = dd;
      for(int h=0; h < constrp(v); h++){
	    d(r) = dd;
	    r++;
      }
    }
	
    // Reduction-step
	d = lambda2 * d;
    alpha = alpha / lambda;
    alpha2 =  alpha*alpha.t();
    lambda = 1;
    lambda2 = 1;
      
    Iteration = iter;
    if (iter > 2){
      if (abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik) {
                                                                              
          break;
        }
      }
    }
    
    
    vec loglik_out;
    int to = Iteration -1;
    loglik_out = loglik.subvec(0, to);
    
    double loglik_max = loglik(to);
    
    
    // if (p * sigma2z < 1e-04){
    //   perror("The estimate of gene expression heritability explained by cis-SNPs is smaller than 0.01%");
    // }
    
    
    
    List output =List::create(Rcpp::Named("alpha") = alpha,
                        Rcpp::Named("sigmaG") = SigGdiag % SigGdiag,
                        Rcpp::Named("sigmaY") = sigma2y,
                        Rcpp::Named("sigmaZdiag") = d,
                        Rcpp::Named("loglik_seq") = loglik_out,
                        Rcpp::Named("loglik") = loglik_max,
                        Rcpp::Named("iteration") = Iteration-1);

   return output;
}// end func


// [[Rcpp::export]]
List GIFT_individualcpppleio(arma::vec X, arma::vec y, arma::mat Zx, arma::mat Zy, arma::mat Omega, arma::vec constralpha, arma::vec constrp, int maxIter, double epsStopLogLik, arma::vec pleioindex){// *
  
  int pxem_indicator = 1;
   	
  // *convert data format
  int k = constrp.n_elem;
  int n2 = y.n_elem, p1 = Zx.n_cols, p2 = Zy.n_cols, n1 = X.n_elem/k;
  if (p1 != p2){
    perror("The dimensions of Zx and Zy are not matched");
  }
  if (k != Omega.n_cols){
    perror("The dimensions of X and Omega are not matched");
  }
  int p = p1;

    // initialization using lmm_pxem
  double sigma2y,loglik0;
  int iter0;
  mat Sigb = zeros<mat>(p,p);
  vec mub  = zeros<vec>(p), d = zeros<vec>(p), SigGdiag = zeros<vec>(k);

  lmm_pxem_ptr2_individual(X, Zx, Omega, maxIter , SigGdiag , d , constrp, k, loglik0 , iter0 , Sigb , mub);

  //declaration of variables used within loop
  mat Omega_inv(k,k), OR(k,k), invOR(k,k), Sigb_inv(p,p), R(p,p), invR(p,p), K1(p,p), K2(p,p), In1 = eye(n1,n1);
  vec mutmp(p), mu(p), res_y(n2), G(n1);
	
  //initialization of parameters
  double lambda = 1, lambda2 = lambda*lambda;
  mat alpha = zeros<mat>(k,1), alpha2(k,k), gamma = zeros<mat>(pleioindex.n_elem,1);
  sigma2y = var(y);

  vec loglik(maxIter), loglik2(maxIter);

  int sum0=0;
  vec constrp1(k+1);
  constrp1(0)= 0;
  for (int v = 0; v < k; v++){
	sum0 += constrp(v);
	constrp1(v+1) = sum0;
  }
	  
  vec contrn=zeros<vec>(k+1);
  int su=0;
  for(int v=0; v < k+1; v++){	  
    contrn(v) = su;
    su += n1;
  }

  OR = chol(Omega);
  invOR = inv(OR);
  Omega_inv = invOR*invOR.t();

  mat K23=zeros<mat>(1,p+1);
  for(int v=0; v < k; v++){	 
	mat K21=zeros<mat>(constrp(v),1);
    if(v!=0){
	  K21 = zeros<mat>(constrp(v),constrp1(v)+1);
    }
    for(int u=v; u < k; u++){	 
      mat K22 = Zy.submat(0,constrp1(v),n2-1,constrp1(v+1)-1).t()*Zy.submat(0,constrp1(u),n2-1,constrp1(u+1)-1);
      K21 = join_rows(K21,K22); 
    }
    K23 =join_cols(K23,K21);
  }
  K23.shed_col(0);
  K23.shed_row(0);
  mat K24 = symmatu(K23);

  mat Ky23=zeros<mat>(1,1);
  for(int v=0; v < k; v++){	 
    mat K22 = Zy.submat(0,constrp1(v),n2-1,constrp1(v+1)-1).t()*y;
    Ky23 = join_cols(Ky23,K22); 
  }
  Ky23.shed_row(0);
	
  mat XK23=zeros<mat>(1,p+1);
  for(int v=0; v < k; v++){	 
    mat XK21=zeros<mat>(constrp(v),1);
    if(v!=0){
	  XK21 = zeros<mat>(constrp(v),constrp1(v)+1);
    }
    for(int u=v; u < k; u++){	 
        mat XK22 = Zx.submat(0,constrp1(v),n1-1,constrp1(v+1)-1).t()*Zx.submat(0,constrp1(u),n1-1,constrp1(u+1)-1) * as_scalar(Omega_inv(v,u));
        XK21 = join_rows(XK21,XK22); 
    }
    XK23 =join_cols(XK23,XK21);
  }
  XK23.shed_col(0);
  XK23.shed_row(0);
  mat XK24 = symmatu(XK23);
  
  mat KX23=zeros<mat>(1,k+1);
  for(int v=0; v < k; v++){	 
  mat KX21=zeros<mat>(constrp(v),1);
  for(int u=0; u < k; u++){	 
  mat KX22 = Zx.submat(0,constrp1(v),n1-1,constrp1(v+1)-1).t()*X.subvec(contrn(u),contrn(u+1)-1);
  KX21 = join_rows(KX21,KX22); 
  }
  KX23 =join_cols(KX23,KX21);
  }
  KX23.shed_col(0);
  KX23.shed_row(0);

  mat Zytu=zeros<mat>(n2,pleioindex.n_elem);
  for(int v=0; v < pleioindex.n_elem; v++){	  
    Zytu.submat(0,v,n2-1,v) = Zy.submat(0,pleioindex(v),n2-1,pleioindex(v));
  }
  gamma = inv(Zytu.t()*Zytu)*Zytu.t()*y;
  
  // initialization of likelihood
  loglik(0) = NAN;
  int Iteration = 1;
  for (int iter = 2; iter <= maxIter; iter ++ ) {
  // E-step

    vec eVal = SigGdiag % SigGdiag;
	mat B1=zeros(1,1);
	vec O2=1.0/SigGdiag;
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=as_scalar(O2(v))*O1;
	  B1=join_cols(B1,O);
    }
    B1.shed_row(0);

    K1 = (B1*B1.t()) % XK24;

	mat A1=zeros(1,1);
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=alpha(v)*O1;
	  A1=join_cols(A1,O);
    }
    A1.shed_row(0);

	K2 = ((A1*A1.t()) % K24)/sigma2y ;

	Sigb_inv = K1 + K2;
    Sigb_inv.diag() += 1.0/d;

    R = chol(Sigb_inv);
    invR = inv(R);
    Sigb = invR*invR.t();

	mat AA = Omega_inv % (O2*O2.t());
	mat tm = KX23 * AA;
    mat in = tm.submat(0,0,constrp(0)-1,0);
    if(k>1){
      for(int u=1; u < k; u++){	 
        mat XK22 = tm.submat(constrp1(u),u,constrp1(u+1)-1,u);
        in = join_cols(in,XK22); 
      }
    }
  
    mat Ky2=zeros<mat>(1,n2);
    for(int v=0; v < k; v++){	 
      mat K2 = Zy.submat(0,constrp1(v),n2-1,constrp1(v+1)-1).t();
      Ky2 = join_cols(Ky2,K2); 
    }
    Ky2.shed_row(0);
    
    mutmp = in + (A1 % (Ky23 - Ky2 * Zytu * gamma))/sigma2y;
   
    mu = Sigb*mutmp;

    double E,M,res_ysum;
  
    res_ysum = n2 -1 - 2*sum(mu % (A1 % Ky23)) + as_scalar(mu.t()*K2*mu)*sigma2y + 2*sum(mu % (A1 %  (Ky2 * Zytu * gamma))) - 2*sum(y % (Zytu * gamma)) + as_scalar((Zytu * gamma).t() * (Zytu * gamma));
  
    E = (n1-1) * accu(Omega % Omega_inv % (O2*O2.t())) - 2*sum(mu % in) + sum(mu%(K1*mu));
	M = accu(mu % (1.0/d)% mu);
      
    loglike_twas_individual(R, E, res_ysum, M, eVal, d, sigma2y, n1 , n2, loglik(iter - 1));
     
    //if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
    //  perror("The likelihood failed to increase!");
    //}
    // M-step
	  
    double tr1 = trace(K1*Sigb);
	   
    if (pxem_indicator == 1){
      lambda = sum(mu % in)/(sum(mu % (K1*mu)) + tr1);
      lambda2 = lambda*lambda;
    }
    else {
      lambda = 1; lambda2 = 1;
    }
    
	mat A2=zeros(n2,1);
    for(int v=0; v < k; v++){	  
      mat O=alpha(v)*Zy.submat(0,constrp1(v),n2-1,constrp1(v+1)-1)*mu.subvec(constrp1(v),constrp1(v+1)-1);
      A2=join_rows(A2,O);
    }
    A2.shed_col(0);
	gamma = inv(Zytu.t()*Zytu)*Zytu.t()*(y - sum(A2,1));
	
    uvec q =  find(constralpha == 0);
	double trde=0;
	double tmp=0;
	mat A = mu*mu.t()+Sigb;
	  
    for(int i = 0;i < k - accu(constralpha);i++){
      int a = q(i);
	  double tmp1=0;
	  for(int j = 0;j < k; j++){
		if(j != a){
          tmp1 += as_scalar(alpha(j)) * trace((K24.submat(constrp1(a),constrp1(j),constrp1(a+1)-1,constrp1(j+1)-1)) * (A.submat(constrp1(j),constrp1(a),constrp1(j+1)-1,constrp1(a+1)-1)));
		}
	  }
    trde = trace((K24.submat(constrp1(a),constrp1(a),constrp1(a+1)-1,constrp1(a+1)-1))*(A.submat(constrp1(a),constrp1(a),constrp1(a+1)-1,constrp1(a+1)-1)));
    tmp = sum(mu.subvec(constrp1(a),constrp1(a+1)-1) % (Ky23.submat(constrp1(a),0,constrp1(a+1)-1,0)));
    alpha(a) = (tmp-tmp1)/trde;
    }
    alpha2 = alpha*alpha.t();

	for(int g = 0; g < k; g++){
	  vec G = X.subvec(contrn(g),contrn(g+1)-1)-lambda*Zx.submat(0,constrp1(g),n1-1,constrp1(g+1)-1)*mu.subvec(constrp1(g),constrp1(g+1)-1);
	  double p1=(sum(G%G)* as_scalar(Omega_inv(g,g)) +lambda2*trace((XK24.submat(constrp1(g),constrp1(g),constrp1(g+1)-1,constrp1(g+1)-1))*(Sigb.submat(constrp1(g),constrp1(g),constrp1(g+1)-1,constrp1(g+1)-1))));
	  double p2=0;
	  for(int u = 0; u < k; u++){
	  vec G1 =  X.subvec(contrn(u),contrn(u+1)-1)-lambda*Zx.submat(0,constrp1(u),n1-1,constrp1(u+1)-1)*mu.subvec(constrp1(u),constrp1(u+1)-1);
	  p2 += ((sum(G%G1)* as_scalar(Omega_inv(g,u)) +lambda2*trace((XK24.submat(constrp1(u),constrp1(g),constrp1(u+1)-1,constrp1(g+1)-1))*(Sigb.submat(constrp1(g),constrp1(u),constrp1(g+1)-1,constrp1(u+1)-1)))))/SigGdiag(u);
	  }
	  SigGdiag(g)= 2*p1/(-(p2-p1/SigGdiag(g))+sqrt(pow((p2-p1/SigGdiag(g)),2)+4*n1*p1));
	}
	
	A1=zeros(1,1);
	for(int v=0; v < k; v++){	  
	  mat O1=ones<mat>(constrp(v),1);
	  mat O=alpha(v)*O1;
	  A1=join_cols(A1,O);
	}
	A1.shed_row(0);
    double tr2 = trace(Sigb*((A1*A1.t()) % K24));
    
    res_ysum = n2 -1 - 2*sum(mu % (A1 % Ky23)) + as_scalar(mu.t()*((A1*A1.t()) % K24)*mu) + 2*sum(mu % (A1 %  (Ky2 * Zytu * gamma))) - 2*sum(y % (Zytu * gamma)) + as_scalar((Zytu * gamma).t() * (Zytu * gamma));
    sigma2y = (res_ysum + tr2)/n2;
    
	vec rr(k);
    int r=0;
	double dd;

	for(int v=0; v < k; v++){
	  dd=sum(mu.subvec(constrp1(v),constrp1(v+1)-1) % mu.subvec(constrp1(v),constrp1(v+1)-1))/constrp(v)+ sum(Sigb.submat(constrp1(v),constrp1(v),constrp1(v+1)-1,constrp1(v+1)-1).diag())/constrp(v);
	  rr(v) = dd;
      for(int h=0; h < constrp(v); h++){
	    d(r) = dd;
	    r++;
      }
    }
	
    // Reduction-step
	d = lambda2 * d;
    alpha = alpha / lambda;
    alpha2 =  alpha*alpha.t();
    lambda = 1;
    lambda2 = 1;
      
    Iteration = iter;
    if (iter > 2){
      if (abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik) {
                                                                              
          break;
        }
      }
    }
    
    
    vec loglik_out;
    int to = Iteration -1;
    loglik_out = loglik.subvec(0, to);
    
    double loglik_max = loglik(to);
    
    
    // if (p * sigma2z < 1e-04){
    //   perror("The estimate of gene expression heritability explained by cis-SNPs is smaller than 0.01%");
    // }
    
    
    
    List output =List::create(Rcpp::Named("alpha") = alpha,
                        Rcpp::Named("sigmaG") = SigGdiag % SigGdiag,
                        Rcpp::Named("sigmaY") = sigma2y,
                        Rcpp::Named("sigmaZdiag") = d,
                        Rcpp::Named("loglik_seq") = loglik_out,
                        Rcpp::Named("loglik") = loglik_max,
                        Rcpp::Named("iteration") = Iteration-1);

   return output;
}// end func
