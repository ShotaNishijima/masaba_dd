
## calculating relative naa (rec=1) at equilibrium ----
calc_rel_naa = function(faa=rep(0,7),M=0.4){
  naa <- 1
  for (i in 1:100) {
    naa = c(naa,as.numeric(rev(naa)[1]*exp(-faa[min(i,length(faa))]-M)))
  }
  tmp = naa[1:(length(faa)-1)]
  c(tmp,sum(naa)-sum(tmp))
}


## predicting weight at age 0
pred_w0 = function(model,Rec,SSB) {
  preddata_w0=model$data %>% 
    mutate(Rec=Rec,SSB=SSB)
  as.numeric(predict(model,newdata=preddata_w0[1,],type="response"))[1]
}

  
## predicting weight at age 1 to A
# model <- mod_w_growth
# w0 = 114.6545
# w0 = opt2$minimum
pred_waa = function(model,w0,naa) {
  weight <- w0
  A <- length(naa)-1
  nc <- naa
  for (j in 1:A) {
    # summary(model_wg)
    preddata_wg = model$data %>% 
      mutate(Age=j,Weight_prev=rev(weight)[1],Cohort_prev=nc[j],Cohort_plus=sum(nc[j:(j+1)]),Number_prev=sum(nc)) %>%
      dplyr::select(Age,Weight_prev,Cohort_prev,Cohort_plus,Number_prev) %>% 
      distinct()
    weight = c(weight,as.numeric(predict(model,newdata=preddata_wg))[1])
  }
  return( weight )
}

delta_mat = function(Maturity,Maturity_prev) (Maturity-Maturity_prev)/(1-Maturity_prev)
y_trans = function(y0,N) (y0*(N-1)+0.5)/N
y0_trans = function(y,N) (y*N-0.5)/(N-1)
trans_mat = function(y0,Maturity_prev) Maturity_prev+y0*(1-Maturity_prev) 

# model <- mod_mat_di
# summary(model)
## predicting maturity at age
pred_mat = function(model,naa,waa,N=nrow(model$model)) {
  nc <- as.numeric(naa)
  maturity <- c(0)
  for(j in 1:3) {
    preddata_mat = model$model %>% 
      mutate(Age=factor(j),Weight=waa[j],Cohort_prev=nc[j],Cohort_plus=sum(nc[j:(j+1)]),Number_prev=sum(nc))
    # tmp = predict(model,newdata=preddata_mat)
    tmp = sort(unique(as.numeric(predict(model,newdata=preddata_mat))))
    tmp = y0_trans(tmp,N)
    tmp = min(1,max(0,tmp))
    tmp = trans_mat(tmp,rev(maturity)[1])
    maturity = c(maturity,tmp)
  }
  maturity <- c(maturity,rep(1,3))
  return( maturity )
}


## calculating weight at age 0 given naa by optimization
opt_w0 = function(naa,model_w0,model_wg,model_mat) {
  nc <- as.numeric(naa)
  # y <- 50
  obj_f2 = function(y,out=FALSE) {
    weight = pred_waa(model_wg,y,nc)
    maturity = pred_mat(model_mat,naa=nc,waa=weight)
    ssb = sum(nc*weight*maturity)
    w0_pred = pred_w0(model_w0,nc[1],ssb)
    if (out==FALSE) {
      return( (y-w0_pred)^2 )
    } else {
      res = data.frame(age=(1:length(nc)-1),waa=weight,maa=maturity,naa=nc,baa=nc*weight,ssb=nc*weight*maturity)
      return( res )
    }
  }
  opt2 = optimize(obj_f2,interval=c(0.1,50))
  # pred_waa(model_wg,opt2$minimum,naa=nc)
  out = obj_f2(opt2$minimum,out=TRUE)
  return( out )
}

## calculating SPR  ----
calc_SPR = function(rec,model_w0,model_wg,model_mat,faa=rep(0,7),M=0.4) {
  rel_naa = calc_rel_naa(faa=faa,M=M)
  naa = rec*rel_naa
  opt = opt_w0(naa,model_w0,model_wg,model_mat)
  opt = opt %>% mutate(faa=faa)
  SSB = sum(opt$ssb)
  Res = list()
  Res$SPRdata = data.frame(R=rec,SSB=SSB,SPR=SSB/rec)
  Res$AGEdata = opt
  return( Res )
}

# resSR=resHS
# 
# resSR=resBH
# resSR$pars
### calculating equilibrium ----
calc_eq = function(faa=rep(0,7),model_w0,model_wg,model_mat,resSR,M=0.4){
  SR = resSR$input$SR
  if (SR == "BH") SRF = function(x,a=resSR$pars$a,b=resSR$pars$b) a*x/(1+b*x)
  if (SR == "HS") SRF = function(x,a=resSR$pars$a,b=resSR$pars$b) min(a*x,a*b)
  obj_fun = function(y) {
    SPRres = calc_SPR(y,model_w0,model_wg,model_mat,faa=faa,M=M)
    ssb = SPRres$SPRdata$SSB
    R = SRF(x=ssb)
    return( (R-y)^2 )
  }
  # y = 5
  # obj_fun(10)
  
  opt = optimize(obj_fun,interval=c(0.000001,50))
  R = opt$minimum
  SPRres = calc_SPR(R,model_w0,model_wg,model_mat,faa=faa,M=M)
  SPR = SPRres$SPRdata$SPR[1]
  A <- nrow(SPRres$AGEdata)-1
  if (sum(faa) == 0) {
    YPR = 0
    SY = 0
  } else {
    refres = frasyr::ref.F(Fcurrent=faa,Pope=TRUE,max.age=100,
                   maa=SPRres$AGEdata$maa,waa=SPRres$AGEdata$waa,M=rep(M,A+1),
                   waa.catch=SPRres$AGEdata$waa,min.age=0,plot=FALSE,pSPR=0,F.range=1)
    YPR = refres$ypr.spr[1,"ypr"]
    SY <- R*YPR
    # if (SR=="BH") SY <- YPR*(resSR$pars$a*SPR-1)/(resSR$pars$b*SPR)
    # if (SR=="HS") {
    #   if (resSR$pars$a*SPR<1) {
    #     SY <- 0
    #   } else {
    #     SY <- YPR*resSR$pars$a*resSR$pars$b
    #   }
    # }
  }
  
  # SPR0res = calc_SPR(R,model_w0,model_wg,model_mat,faa=faa*0,M=M)
  # refres$ypr.spr
  # SPRres$AGEdata
  res1 = SPRres$SPRdata %>%
    mutate(YPR = YPR, B = sum(SPRres$AGEdata$baa), SY=SY, 
           SPR0 =   refres$spr0) %>%
    mutate(pSPR = SPR/SPR0)
  res2 = SPRres$AGEdata %>% 
    mutate(faa=faa,M=M) %>%
    mutate(catch=baa*exp(-0.5*M)*(1-exp(-faa)))
  return( list(EQdata=res1,AGEdata=res2))
}


calc_eq_fix = function(faa=rep(0,7),waa,maa,resSR,M=0.4){
  SR = resSR$input$SR
  if (SR == "BH") SRF = function(x,a=resSR$pars$a,b=resSR$pars$b) a*x/(1+b*x)
  if (SR == "HS") SRF = function(x,a=resSR$pars$a,b=resSR$pars$b) min(a*x,a*b)
  A <- length(faa)-1
  obj_fun = function(y,out=FALSE) {
    naa = y*calc_rel_naa(faa,M=M)
    baa = naa*waa
    ssb = baa*maa
    SSB = sum(ssb)
    R = SRF(x=SSB)
    if (out==FALSE) return( (R-y)^2 )
    if (out==TRUE) {
      res1 = data.frame(age=0:A,waa=waa,maa=maa,naa=naa,baa=baa,ssb=ssb)
      res2 = data.frame(R=R,SSB=SSB,SPR=SSB/R)
      return( list(AGEdata=res1,SPRdata=res2) )
    }
  }
  opt = optimize(obj_fun,interval=c(0.000001,50))
  R = opt$minimum
  eqres = obj_fun(R,out=TRUE)
  if (sum(faa) == 0) {
    YPR = 0
    SY = 0
  } else {
    refres = frasyr::ref.F(Fcurrent=faa,Pope=TRUE,max.age=100,
                           maa=maa,waa=waa,M=rep(M,A+1),
                           waa.catch=waa,min.age=0,plot=FALSE,pSPR=0,F.range=1)
    YPR = refres$ypr.spr[1,"ypr"]
    SY <- R*YPR
  }
  
  res1 = eqres$SPRdata %>%
    mutate(YPR = YPR, B = sum(eqres$AGEdata$baa), SY=SY, 
           SPR0 =   refres$spr0) %>%
    mutate(pSPR = SPR/SPR0)
  res2 = eqres$AGEdata %>% 
    mutate(faa=faa,M=M) %>%
    mutate(catch=baa*exp(-0.5*M)*(1-exp(-faa)))
  return( list(EQdata=res1,AGEdata=res2))
}
