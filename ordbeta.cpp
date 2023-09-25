#include <TMB.hpp>
#include <iostream>

template <class Type>
Type square(Type x){return x*x;}

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
    ans = -logspace_add(Type(0), -eta);  //log(1/(1+exp(-eta))) = -log(1+exp(-eta))
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
  DATA_ARRAY(obs_g);
  PARAMETER_VECTOR(alpha);
  PARAMETER(logdisp);
  // PARAMETER_VECTOR(psi);
  PARAMETER(omicron);
  PARAMETER_ARRAY(logitmaa);
  
  // vector<Type> y_pred = invlogit(alpha(x));
  Type ans=0.0;
  Type disp=exp(logdisp);
  Type s1,s2;
  // Type s3;
  Type v,r;
  array<Type> maa_true(logitmaa.rows(),logitmaa.cols());
  array<Type> maa_pred(logitmaa.rows(),logitmaa.cols());
  
  // process model
  for(int j=0;j<logitmaa.cols();j++){
    for(int i=0;i<logitmaa.rows();i++){
      maa_true(i,j) = invlogit(logitmaa(i,j));
    }
  }
  
  for(int j=1;j<logitmaa.cols();j++){
    for(int i=0;i<logitmaa.rows();i++){
      if(i==0){
        maa_pred(i,j)=invlogit(alpha(i));
      }else{
        maa_pred(i,j)=maa_true(i-1,j-1)+invlogit(alpha(i))*(Type(1.0)-maa_true(i-1,j-1));
      }
      // ans+=-dnorm(logit(maa_true(i,j)),logit(maa_pred(i,j)),exp(omicron),true);
      s1=maa_pred(i,j)*exp(omicron);
      s2=(Type(1.0)-maa_pred(i,j))*exp(omicron);
      s1=maa_pred(i,j)*disp; //assuming observation error size = measurement error size
      s2=(Type(1.0)-maa_pred(i,j))*disp;
      v=pbeta(maa_true(i,j),s1,s2);
      r=qnorm(v,Type(0.0),Type(1.0));
      ans+=-dnorm(r,Type(0.0),Type(1.0));
    }
  }

  // observation model
  int a,y;
  for(int i=0;i<obs_g.rows();i++){
    a=CppAD::Integer(obs_g(i,0)-obs_g(0,0));
    y=CppAD::Integer(obs_g(i,1)-obs_g(0,1));
    // if (obs_g(i,2) == 0.0) {
    //   ans += -log1m_inverse_linkfun(logit(maa_true(a,y)) - psi(0), logit_link);
    //   // std::cout << "zero " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(tmp_loglik) << std::endl;
    // } else if (obs_g(i,2) == 1.0) {
    //   ans += -log_inverse_linkfun(logit(maa_true(a,y)) - psi(1), logit_link);
    //   // std::cout << "one " << asDouble(eta(i)) << " " << asDouble(psi(1)) << " " << asDouble(tmp_loglik) << std::endl;
    // } else {
      s1 = maa_true(a,y)*disp;
      s2 = (Type(1)-maa_true(a,y))*disp;
      // s3 = logspace_sub(log_inverse_linkfun(logit(maa_true(a,y)) - psi(0), logit_link),
      //                   log_inverse_linkfun(logit(maa_true(a,y)) - psi(1), logit_link));
      // ans += -s3 - dbeta(obs_g(i,2), s1, s2, true);
      ans += -dbeta(obs_g(i,2), s1, s2, true);
    // }

    // s1 = y_pred(i)*disp;
    // s2 = (Type(1)-y_pred(i))*disp;
    // ans += -dbeta(y(i), s1, s2, true);
    // }
    // ans+=-dnorm(obs_g(i,2),maa_true(a,y),disp,true);
  }
  
  // ADREPORT(alpha);
  
  return ans;
}