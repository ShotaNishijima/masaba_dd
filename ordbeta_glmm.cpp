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
  DATA_VECTOR(y);
  DATA_IVECTOR(id);
  
  PARAMETER(logitp);
  PARAMETER_VECTOR(logitr);
  PARAMETER_VECTOR(logdisp);
  PARAMETER_VECTOR(psi);
  
  Type ans=0.0;
  vector<Type> disp=exp(logdisp);
  vector<Type> y_pred(y.size());
  Type s1,s2,s3;
  // Type v,w;
  vector<Type> r=invlogit(logitr);
  
  for(int i=0;i<r.size();i++){
    // s1=invlogit(logitp)*disp(0);
    // s2=(Type(1.0)-invlogit(logitp))*disp(0);
    // v=pbeta(r(i),s1,s2);
    // w=qnorm(v,Type(0.0),Type(1.0));
    // ans += -dnorm(w,Type(0.0),Type(1.0),true); // random effect component
    ans += -dnorm(logitr(i),logitp,disp(0),true); // random effect component
  }
  
  for(int i=0;i<y.size();i++){
    // y_pred(i) = logitp + r(id(i));
    // y_pred(i) = invlogit(y_pred(i));
    // ans += -dnorm(y(i),y_pred(i),disp(1),true);
    y_pred(i)=r(id(i));
    if (y(i) == 0.0) {
      ans += -log1m_inverse_linkfun(logit(y_pred(i)) - psi(0), logit_link);
      // std::cout << "zero " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(tmp_loglik) << std::endl;
    } else if (y(i) == 1.0) {
      ans += -log_inverse_linkfun(logit(y_pred(i)) - psi(1), logit_link);
      // std::cout << "one " << asDouble(eta(i)) << " " << asDouble(psi(1)) << " " << asDouble(tmp_loglik) << std::endl;
    } else {
      s1 = y_pred(i)*disp(1);
      s2 = (Type(1.0)-y_pred(i))*disp(1);
      s3 = logspace_sub(log_inverse_linkfun(logit(y_pred(i)) - psi(0), logit_link),
                        log_inverse_linkfun(logit(y_pred(i)) - psi(1), logit_link));
      ans += -s3 - dbeta(y(i), s1, s2, true);
    }
  }
  
  return ans;
  
}