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
  DATA_ARRAY(obs_g); // observation matrix of maturity (column 0: age, 1: year, 2: maturity (0-1))
  DATA_ARRAY(naa);
  // DATA_ARRAY(maa);

  PARAMETER_ARRAY(logwaa);
  PARAMETER_ARRAY(logitmaa);
  // PARAMETER_ARRAY(r); // latent variables from standard normal distribution for maturity process

  PARAMETER_VECTOR(beta_w0);
  PARAMETER_VECTOR(alpha_w);
  PARAMETER_VECTOR(rho_w);
  PARAMETER(iota); // log(SD) for process error in weight growth
  PARAMETER(omicron); // log(SD) for process error in maturity growth
  PARAMETER(logCV_w); // CV for observation error in weight
  PARAMETER_VECTOR(alpha_g); // intercept of maturity modeling
  PARAMETER_VECTOR(psi);
  PARAMETER(logdisp); // dispersion parameter (phi) in log space

  // int timeSteps=waa.dim[1]; // 年数
  // int stateDimN=waa.dim[0]; // number of age classes

  DATA_INTEGER(minAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_INTEGER(dist_wobs); // probability distribution for observation of weight (0: lognormal, 1: gamma)
  DATA_SCALAR(scale_number);
  DATA_VECTOR(g_fix); // its length is the number of age classes. Non-negative value (0-1) represents the fixed value of maturity at age while a negative value (e.g., -1)indicates estimating maturity.

  DATA_SCALAR(scale);
  DATA_VECTOR(ssb);

  Type sd_w=exp(iota);
  // Type sd_g=exp(omicron);

  array<Type> logwaa_pred(logwaa.rows(),logwaa.cols());
  array<Type> waa_true(logwaa.rows(),logwaa.cols());
  array<Type> maa_true(logwaa.rows(),logwaa.cols());
  array<Type> maa_pred(logitmaa.rows(),logitmaa.cols());

  Type ans_w=0;
  vector<Type> wp(2);
  Type alpha_w_total, rho_w_total,beta_w0_total;
  vector<Type> N_sum(logwaa.cols());

  // process model for weight
  for(int j=0;j<logwaa.cols();j++){
    N_sum(j)=Type(0.0);
    for(int i=0;i<logwaa.rows();i++){
      waa_true(i,j)=exp(logwaa(i,j));
      N_sum(j)+=naa(i,j)/scale_number;
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
  Type multi_g;
  Type ans_g=0.0;

  for(int j=0;j<logwaa.cols();j++){
    // N_sum(j)=Type(0.0);
    for(int i=0;i<logwaa.rows();i++){
      if(g_fix(i)>-0.5){
        maa_true(i,j)=g_fix(i);
      }else{
        a=CppAD::Integer(-g_fix(i))-1;
        maa_true(i,j)=invlogit(logitmaa(a,j));
        // g_pred(a,j)
      }
      // N_sum(j)+=naa(i,j)/scale_number;
    }
  }

  // process likelihood
  Type v,r;
  Type s1, s2;
  for(int j=1;j<logwaa.cols();j++){ //最初の年は除く（2年目から）
    for(int i=0;i<logwaa.rows();i++){
      if(g_fix(i)<0.0){
        a=CppAD::Integer(-g_fix(i))-1;
        multi_g=alpha_g(a);
        maa_pred(a,j)=maa_true(i-1,j-1)+invlogit(multi_g)*(Type(1.0)-maa_true(i-1,j-1));
        s1=maa_pred(a,j)*exp(omicron);
        s2=(Type(1.0)-maa_pred(a,j))*exp(omicron);
        v=pbeta(maa_true(i,j),s1,s2);
        r=qnorm(v,Type(0.0),Type(1.0));
        ans_g += -dnorm(r,Type(0.0),Type(1.0),true);  
        // logitg_pred(a,j)=invlogit(multi_g); // increasing ratio of maturity (未成熟だった個体のうち成熟した割合)
        // logitg_pred(a,j)*=(Type(1.0)-maa_true(i-1,j-1));
        // logitg_pred(a,j)+=maa_true(i-1,j-1); //expected maturity
        // logitg_pred(a,j)=logit(logitg_pred(a,j));
        // ans_g+=-dnorm(logitg(a,j),logitg_pred(a,j),sd_g,true); // assuming logit-normal process error
      }
    }
  }
   
  // observation model for maturity
  // int link=1; //logit link
  minYear=CppAD::Integer((obs_g(0,1)));
  int minAge_g=CppAD::Integer((obs_g(0,0))); 
  Type s3;
  // Type tmp_loglik;
  vector<Type> eta(obs_g.rows());
  Type disp=exp(logdisp);

  for(int i=0;i<obs_g.rows();i++){
    a=CppAD::Integer(obs_g(i,0))-minAge_g;
    y=CppAD::Integer(obs_g(i,1))-minYear;
    
    if (obs_g(i,2) == 0.0) {
      // ans_g += -log1m_inverse_linkfun(alpha_g(a) - psi(0), logit_link);
      ans_g += -log1m_inverse_linkfun(logit(maa_true(a,y)) - psi(0), logit_link);
      // std::cout << "zero " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(tmp_loglik) << std::endl;
    } else if (obs_g(i,2) == 1.0) {
      // ans_g += -log_inverse_linkfun(alpha_g(a) - psi(1), logit_link);
      ans_g += -log_inverse_linkfun(logit(maa_true(a,y)) - psi(1), logit_link);
      // std::cout << "one " << asDouble(eta(i)) << " " << asDouble(psi(1)) << " " << asDouble(tmp_loglik) << std::endl;
    } else {
      // s1 = invlogit(alpha_g(a))*disp;
      s1 = maa_true(a,y)*disp;
      s2 = (Type(1)-maa_true(a,y))*disp;
      // s3 = logspace_sub(log_inverse_linkfun(alpha_g(a) - psi(0), logit_link),
      //                   log_inverse_linkfun(alpha_g(a) - psi(1), logit_link));
      s3 = logspace_sub(log_inverse_linkfun(logit(maa_true(a,y)) - psi(0), logit_link),
                        log_inverse_linkfun(logit(maa_true(a,y)) - psi(1), logit_link));
      ans_g += -s3 - dbeta(obs_g(i,2), s1, s2, true);
    }
  }

  // for(int i=0;i<obs_g.rows();i++){
  //   a=CppAD::Integer(obs_g(i,0))-minAge_g;
  //   y=CppAD::Integer(obs_g(i,1))-minYear;
  //   // eta(i)=logit(maa_true(a,y));
  //   eta(i)=alpha_g(a);
  //   
  //   //use ordbeta family in glmmTMB
  //   // https://github.com/glmmTMB/glmmTMB/blob/a74c35fa443c9e86698676f3ffc31149ab1df849/glmmTMB/src/glmmTMB.cpp#L676C2-L689C1
  //   // https://github.com/saudiwin/ordbetareg_pack/blob/master/R/modeling.R#L565-L573
  //   
  //   if(obs_g(i,2) == 0.0){
  //     ans_g += -log1m_inverse_linkfun(eta(i) - psi(0), log_link);
  //     // std::cout << "zero " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(tmp_loglik) << std::endl;
  //   }else
  //     if(obs_g(i,2) == 1.0){
  //     ans_g += -log_inverse_linkfun(eta(i) - psi(1), log_link);
  //     // std::cout << "one " << asDouble(eta(i)) << " " << asDouble(psi(1)) << " " << asDouble(tmp_loglik) << std::endl;
  //   }else{
  //     // s1 = maa_true(a,y)*disp;
  //     // s2 = (Type(1)-maa_true(a,y))*disp;
  //     s1 = invlogit(alpha_g(a))*disp;
  //     s2 = (Type(1)-invlogit(alpha_g(a)))*disp;
  //     s3 = logspace_sub(log_inverse_linkfun(eta(i) - psi(0), log_link),
  //                       log_inverse_linkfun(eta(i) - psi(1), log_link));
  //     ans_g += -s3 - dbeta(obs_g(i,2), s1, s2, true);
  //     }
  //   // ans_g+=-tmp_loglik;
  // }

  Type ans=ans_w+ans_g;

  ADREPORT(waa_true);
  ADREPORT(maa_true);

  REPORT(logwaa);
  REPORT(sd_w);
  // REPORT(sd_g);
  REPORT(shape);
  REPORT(logwaa_pred);
  // REPORT(logitg_pred);
  REPORT(disp);
  REPORT(eta);
  REPORT(ans_w);
  REPORT(ans_g);

  return ans;
}
