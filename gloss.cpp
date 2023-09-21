// State space assessment model of growth -------------------------------------

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

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(obs_w); // observation matrix (column 0: age, 1: year, 2: weight (g))
  // DATA_IVECTOR(ss_wg); // whether state-space modeling is applied for weight and maturity
  // DATA_ARRAY(obs_g);
  DATA_ARRAY(naa);
  // DATA_ARRAY(maa);

  PARAMETER_ARRAY(logwaa);
  // PARAMETER_ARRAY(maa);

  PARAMETER_VECTOR(beta_w0);
  PARAMETER_VECTOR(alpha_w);
  PARAMETER_VECTOR(rho_w);
  PARAMETER(iota); // log(SD) for process error in weight growth
  // PARAMETER(omicron); // log(SD) for process error in maturity growth
  PARAMETER(logCV_w); // CV for observation error in weight


  // int timeSteps=waa.dim[1]; // 年数
  // int stateDimN=waa.dim[0]; // number of age classes

  DATA_INTEGER(minAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_INTEGER(dist_wobs); // probability distribution for observation of weight (0: lognormal, 1: gamma)

  Type sd_w = exp(iota);

  array<Type> logwaa_pred(logwaa.rows(),logwaa.cols());
  array<Type> waa_true(logwaa.rows(),logwaa.cols());
  Type ans_w=0;
  vector<Type> wp(2);
  Type alpha_w_total, rho_w_total,beta_w0_total;

  // process model for weight
  for(int i=0;i<logwaa.rows();i++){
    for(int j=0;j<logwaa.cols();j++){
      waa_true(i,j)=exp(logwaa(i,j));
    }
  }

  for(int j=1;j<logwaa.cols();j++){ //最初の年は除く（2年目から）
    for(int i=0;i<logwaa.rows();i++){
      alpha_w_total=alpha_w(0);
      rho_w_total=rho_w(0);
      beta_w0_total=beta_w0(0);
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

  ADREPORT(waa_true);

  REPORT(logwaa);
  REPORT(sd_w);
  REPORT(shape);

  return ans_w;
}
