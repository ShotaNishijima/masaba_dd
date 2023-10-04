// State space assessment model with growth -------------------------------------

#include <TMB.hpp>
#include <iostream>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type>
Type square(Type x){return x*x;}

// sqrt
template<class Type>
Type sqrt(Type x){
  return pow(x,Type(0.5));
}

// inverse of logit
// template<class Type>
// Type invlogit(Type x){
//   return Type(1.0)/(Type(1.0)+exp(-x));
// }

// https://github.com/glmmTMB/glmmTMB/blob/a74c35fa443c9e86698676f3ffc31149ab1df849/glmmTMB/src/glmmTMB.cpp#L50C1-L58C3
enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  cloglog_link             = 4,
  identity_link            = 5,
  sqrt_link                = 6
};

// https://github.com/glmmTMB/glmmTMB/blob/a74c35fa443c9e86698676f3ffc31149ab1df849/glmmTMB/src/glmmTMB.cpp#L81C1-L111C2
template<class Type>
Type inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = exp(eta);
    break;
  case identity_link:
    ans = eta;
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case cloglog_link:
    ans = Type(1) - exp(-exp(eta));
    break;
  case inverse_link:
    ans = Type(1) / eta;
    break;
  case sqrt_link:
    ans = eta*eta; // pow(eta, Type(2)) doesn't work ... ?
    break;
    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

// https://github.com/glmmTMB/glmmTMB/blob/a74c35fa443c9e86698676f3ffc31149ab1df849/glmmTMB/src/glmmTMB.cpp#L134C1-L149C2
/* log transformed inverse_linkfun without losing too much accuracy */
template<class Type>
Type log_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = eta;
    break;
  case logit_link:
    ans = -logspace_add(Type(0), -eta); //log(1/(1+exp(-eta))) = -log(1+exp(-eta))
    break;
  default:
    ans = log( inverse_linkfun(eta, link) );
  } // End switch
  return ans;
}


// https://github.com/glmmTMB/glmmTMB/blob/a74c35fa443c9e86698676f3ffc31149ab1df849/glmmTMB/src/glmmTMB.cpp#L151C1-L166C2
/* log transformed inverse_linkfun without losing too much accuracy */
template<class Type>
Type log1m_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = logspace_sub(Type(0), eta);
    break;
  case logit_link:
    ans = -logspace_add(Type(0), eta);
    break;
  default:
    ans = logspace_sub(Type(0), log( inverse_linkfun(eta, link) ));
  } // End switch
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(obs_w); // observation matrix of weight (column 0: age, 1: year, 2: weight (g))
  // DATA_IVECTOR(ss_wg); // whether state-space modeling is applied for weight and maturity
  DATA_ARRAY(maa);
  DATA_ARRAY(naa);

  PARAMETER_ARRAY(logwaa);
  // PARAMETER_ARRAY(logitmaa);
  PARAMETER_ARRAY(r); // latent variables from standard normal distribution for maturity process

  PARAMETER_VECTOR(beta_w0);
  PARAMETER_VECTOR(alpha_w);
  PARAMETER_VECTOR(rho_w);
  PARAMETER(omicron); // log(SD) for process error in weight growth
  PARAMETER(iota); // log(SD) for process error in maturity growth
  PARAMETER(logCV_w); // CV for observation error in weight
  PARAMETER_VECTOR(alpha_g); // intercept of maturity modeling
  PARAMETER_VECTOR(psi);
  PARAMETER(logdisp); // dispersion parameter (phi) in log space for beta distribution

  // int timeSteps=waa.dim[1]; // 年数
  // int stateDimN=waa.dim[0]; // number of age classes

  DATA_INTEGER(minAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_INTEGER(dist_wobs); // probability distribution for observation of weight (0: lognormal, 1: gamma)
  DATA_SCALAR(scale_number);
  DATA_VECTOR(g_fix); // its length is the number of age classes. Non-negative value (0-1) represents the fixed value of maturity at age while a negative value (e.g., -1)indicates estimating maturity.

  DATA_SCALAR(scale);
  DATA_VECTOR(ssb);
  // vector<Type> ssb(naa.cols());

  Type sd_w=exp(omicron);
  Type sd_g=exp(iota);

  array<Type> logwaa_pred(logwaa.rows(),logwaa.cols());
  array<Type> waa_true(logwaa.rows(),logwaa.cols());
  logwaa_pred.fill(0.0);
  waa_true.fill(-1);

  Type ans_w=0;
  vector<Type> wp(2);
  Type alpha_w_total, rho_w_total,beta_w0_total;
  vector<Type> N_sum(logwaa.cols());

  // process model for weight
  for(int j=0;j<logwaa.cols();j++){
    N_sum(j)=Type(0.0);
    // ssb(j)=Type(0.0);
    for(int i=0;i<logwaa.rows();i++){
      waa_true(i,j)=exp(logwaa(i,j));
      N_sum(j)+=naa(i,j)/scale_number;
      // ssb(j)+=naa(i,j)*maa(i,j)*waa_true(i,j)/scale;
    }
  }

  for(int j=1;j<logwaa.cols();j++){ //最初の年は除く（2年目から）
    for(int i=0;i<logwaa.rows();i++){
      alpha_w_total=alpha_w(0);
      if(alpha_w.size()>1){
        alpha_w_total+=alpha_w(1)*N_sum(j);
      }
      rho_w_total=rho_w(0);
      beta_w0_total=beta_w0(0);
      if(beta_w0.size()>1){
        beta_w0_total+=beta_w0(1)*ssb(j)/scale;
      }
      if(i==0){ // age 0
        logwaa_pred(i,j)=beta_w0_total;
      }else{
        if(i<logwaa.rows()-1){
          // from Brody-growth coefficient and the von-Bertalanffy weight model
          logwaa_pred(i,j)=alpha_w_total;
          logwaa_pred(i,j)+=rho_w_total*waa_true(i-1,j-1);
          logwaa_pred(i,j)=log(logwaa_pred(i,j));
        }else{
          if(maxAgePlusGroup==1) { //plus group
            wp(0)=alpha_w_total;
            wp(0)+=rho_w_total*waa_true(i-1,j-1);
            wp(1)=alpha_w_total;
            wp(1)+=rho_w_total*waa_true(i,j-1);
            logwaa_pred(i,j)=naa(i-1,j-1)*wp(0)+naa(i,j-1)*wp(1);
            logwaa_pred(i,j)=logwaa_pred(i,j)/(naa(i-1,j-1)+naa(i,j-1));
            logwaa_pred(i,j)=log(logwaa_pred(i,j));
          }else{
            logwaa_pred(i,j)=alpha_w_total;
            logwaa_pred(i,j)+=rho_w_total*waa_true(i-1,j-1);
            logwaa_pred(i,j)=log(logwaa_pred(i,j));
          }
        }
      }
      ans_w+=-dnorm(logwaa(i,j),logwaa_pred(i,j),sd_w,true);
    }
  }

  // observation model for weight
  int minYear=CppAD::Integer((obs_w(0,1)));
  int a, y;
  Type scale_par;
  Type shape=pow(exp(logCV_w),Type(-2.0));
  if (dist_wobs==0) { //lognormal
    shape=sqrt(log(Type(1.0)+pow(exp(logCV_w),Type(2.0)))); //SD for lognormal distribution
  }
  for(int i=0;i<obs_w.rows();i++){
    a=CppAD::Integer(obs_w(i,0))-minAge;
    y=CppAD::Integer(obs_w(i,1))-minYear;
    scale_par=waa_true(a,y)/shape;
    if (dist_wobs==1){
      ans_w+=-dgamma(obs_w(i,2),shape,scale_par,true); //using Gamma distribution
    }else{ // lognormal
      ans_w+=-dnorm(log(obs_w(i,2)),logwaa(a,y)-Type(0.5)*shape*shape,shape,true); //using lognormal distribution with bias correction
    }
  }

  // process model for maturity (no observation error)
  Type ans_g=0.0;
  array<Type> maa_pred(maa.rows(),maa.cols());
  array<Type> maa_diff(maa.rows(),maa.cols());
  array<Type> maa_true(maa.rows(),maa.cols());
  array<Type> diff_pred(maa.rows(),maa.cols());
  maa_pred.fill(-1);
  maa_diff.fill(-1);
  Type multi_g;

  Type disp=exp(logdisp);
  Type s1, s2, s3;
  for(int j=1;j<maa.cols();j++){ //最初の年は除く（2年目から）
    for(int i=1;i<maa.rows();i++){
      maa_diff(i,j)=(maa(i,j)-maa(i-1,j-1))/(Type(1.0)-maa(i-1,j-1));
      if(g_fix(i)<0.0){
        a=CppAD::Integer(-g_fix(i))-1;
        multi_g=alpha_g(a);
        maa_pred(i,j)=maa(i-1,j-1)+invlogit(multi_g)*(Type(1.0)-maa(i-1,j-1));
        if (maa_diff(i,j) == 0.0) {
          // ans_g += -log1m_inverse_linkfun(logit(maa_pred(i,j)) - psi(0), logit_link);
          ans_g += -log1m_inverse_linkfun(multi_g - psi(0), logit_link);
          // std::cout << "zero " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(tmp_loglik) << std::endl;
        } else if (maa_diff(i,j) == 1.0) {
          // ans_g += -log_inverse_linkfun(logit(maa_pred(i,j)) - psi(1), logit_link);
          ans_g += -log_inverse_linkfun(multi_g - psi(1), logit_link);
          // std::cout << "one " << asDouble(eta(i)) << " " << asDouble(psi(1)) << " " << asDouble(tmp_loglik) << std::endl;
        } else {
          // s1 = maa_pred(i,j)*disp;
          // s2 = (Type(1)-maa_pred(i,j))*disp;
          s1 = invlogit(multi_g)*disp;
          s2 = (Type(1.0)-invlogit(multi_g))*disp;
          // s3 = logspace_sub(log_inverse_linkfun(logit(maa_pred(i,j)) - psi(0), logit_link),
          //                   log_inverse_linkfun(logit(maa_pred(i,j)) - psi(1), logit_link));
          s3 = logspace_sub(log_inverse_linkfun(multi_g - psi(0), logit_link),
                            log_inverse_linkfun(multi_g - psi(1), logit_link));
          ans_g += -s3 - dbeta(maa_diff(i,j), s1, s2, true);
        }
      }
    }
  }
   
  Type ans=ans_w+ans_g;

  ADREPORT(waa_true);
  // ADREPORT(maa_true);

  REPORT(logwaa);
  REPORT(sd_w);
  REPORT(shape);
  REPORT(logwaa_pred);
  // REPORT(maa_pred);
  // REPORT(disp);
  REPORT(ans_w);
  // REPORT(ans_g);

  return ans;
}
