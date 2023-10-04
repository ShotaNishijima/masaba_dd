

### set working directory ---- 

gitdir = "~/git/masaba_dd"
setwd(gitdir)

savedir = "res"

if (!file.exists(savedir)) dir.create(savedir)

savename = function(x) paste0(savedir,"/",x)

### read packages ----

# devtools::install_github("ichimomo/frasyr", ref="dev")

library(frasyr)
library(tidyverse)
library(MASS)
library(MuMIn)
options(na.action = "na.fail")
library("dplyr", character.only = TRUE)

# install.packages("broom.mixed")
# install.packages("betareg")
library(broom.mixed)
library(betareg)
# library(lmtest)
# install.packages("effects")
# library(effects)
# install.packages("investr")
# library(investr)
library(glmmTMB)
library(TMB)

getwd()


### glmm via TMB with ordbeta family

logit = function(x) log(x/(1-x))
# logit = Vectorize(logit)

logit(c(0.1,0.2))

invlogit = function(x) 1/(1+exp(-x))

id = rep(1:10,10)
# disp = c(2,2)
disp = c(0.1,0.1)
p = 0.3
set.seed(1)
r = rbeta(10,p*disp[1],(1-p)*disp[1])
hist(r)
y0 = sapply(1:100, function(i) rbeta(1,r[id[i]]*disp[2],(1-r[id[i]])*disp[2]))
y=round(y0,2)
y
hist(y)

data = data.frame(y=y,id=factor(id))

res =glmmTMB(y~1+(1|id),data=data,family=ordbeta)
summary(res)

predict(res) %>% unique()

# r = ranef(res)$cond$id

res$sdr$par.random %>% mean

mean(unlist(r))

logit(p)

compile("ordbeta_glmm.cpp")
dyn.load(dynlib("ordbeta_glmm"))

tmbdata = list(y=y,id=id-1)

# boxplot(y~id)

pars = list(logitp=logit(p),logitr=rep(0,10),logdisp=log(disp),
            psi=as.numeric(res$fit$par[names(res$fit$par)=="psi"]))

# pars = list(logitp=logit(p),r=rep(0,10),logdisp=log(disp))
# pars = list(logitp=logit(p),r=r-mean(r),logdisp=c(0,0.1))

pars$logitp %>% class

obj = TMB::MakeADFun(tmbdata,pars,random=c("logitr"),DLL="ordbeta_glmm")


# map = list()

# pars$maa
# 
# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss",
#                      map=map)
# 
# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss")

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

rep = sdreport(obj,bias.correct = TRUE)
rep



## 正規分布のランダム効果だと上手く推定できる

# 
# compile("ordbeta.cpp")
# dyn.load(dynlib("ordbeta"))
# 
# 
# tmbdata = list(obs_g=obs_g)
# 
# alpha_g = mod_mat0$fit$par[names(mod_mat0$fit$par)=="beta"] %>% as.numeric()
# psi = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="psi"])
# logdisp = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="betad"])[1]
# 
# maa = vpares$input$dat$maa
# 
# logitmaa = matrix(
#   rep(logit(rowMeans(maa[2:4,])),ncol(logwaa)),nrow=3)
# 
# naa0 = vpares$naa/1000
# 
# naa = naa0[-1,]+naa0[-nrow(naa0),] %>% as.matrix()
# 
# naa=as.matrix(naa[1:3,])
# # cohort_plus
# 
# class(obs_g)
# class(naa)
# 
# tmbdata = list(obs_g=obs_g,naa=naa)
# 
# # as.matrix(tmbdata$naa)
# 
# p=0.3
# disp=10000
# hist(rbeta(1000,p*disp,(1-p)*disp))
# 
# pars = list(alpha=alpha_g,logdisp=logdisp,psi=psi,omicron=log(1000),
#             logitmaa=logitmaa,beta=rep(0,length(alpha_g)))
# 
# pars = list(alpha=alpha_g,logdisp=log(1000),psi=psi,omicron=1,
#             logitmaa=logitmaa,beta=rep(0,length(alpha_g)))
# 
# 
# map=list()
# # map = list(omicron=factor(NA))
# map = list(logdisp=factor(NA))
# 
# # map = list(omicron=factor(NA),psi=rep(factor(NA),2))
# 
# obj = TMB::MakeADFun(tmbdata,pars,random=c("logitmaa"),DLL="ordbeta",
#                      map=map)
# 
# # obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss_orig")
# 
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# opt
# 
# rep = sdreport(obj,bias.correct = TRUE)
# rep
# 
# 
# 


vpares = get(load("data/vpa_masaba_P2022.rda")) # assessment result in 2022

range(colSums(vpares$ssb))
range(colSums(vpares$ssb))/1000

waa_obs = vpares$input$dat$waa 

waa_dat = expand.grid(Age=as.numeric(rownames(vpares$naa)),
                      Year=as.numeric(colnames(vpares$naa))) %>% 
  mutate(Weight=as.numeric(unlist(vpares$input$dat$waa)),
         Maturity=as.numeric(unlist(vpares$input$dat$maa)))

A = waa_dat$Age %>% max

waa_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if (waa_dat$Age[i]<A) {
      value <- as.numeric(
        vpares$input$dat$waa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
    } else {
      value1 <- as.numeric(
        vpares$input$dat$waa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      value2 <- as.numeric(
        vpares$input$dat$waa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)]) # plus group
      n1 <- as.numeric(
        vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      n2 <- as.numeric(
        vpares$naa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)]) # plus group
      value <- (value1*n1+value2*n2)/(n1+n2)
    }
  }
  value
})

maa_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    value <- as.numeric(
      vpares$input$dat$maa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
  }
  value
})

# Total number (from age 0 to 6+)
num_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    value <- as.numeric(
      sum(vpares$naa[,as.character(waa_dat$Year[i]-1)]))
  }
  value
})

# picking up the same cohort(同じコホートにおける前年の資源尾数)
cohort_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if(waa_dat$Age[i]==0) {
      value <- NA
    } else {
      # if (waa_dat$Age[i]<A) {
      value <- as.numeric(vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      # } else {
      # }
    }
  }
  value
})

# 前年の同じコホートと1歳上のコホートの資源尾数の合計
cohort_plus = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if(waa_dat$Age[i]==0) {
      value <- NA
    } else {
      # if (waa_dat$Age[i]<A) {
      value <- as.numeric(vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      value <- value + as.numeric(vpares$naa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)])
      # } else {
      # }
    }
  }
  value
})

# scale_number = 1000

waa_dat = waa_dat %>%
  mutate(Weight_prev=waa_prev,
         Maturity_prev=maa_prev,
         Number_prev=num_prev/1000,
         Cohort_prev=cohort_prev/1000,
         Cohort_plus=cohort_plus/1000) %>%
  mutate(logN_prev=log(num_prev)) %>%
  mutate(YearClass=Year-Age) %>%
  # mutate(interact_norm = Weight_prev*Number_prev,
  #        interact_log = Weight_prev*logN_prev) %>%
  filter(Year<max(Year))


# waa_dat = expand.grid(Age=as.numeric(rownames(vpares$naa)),
#                       Year=as.numeric(colnames(vpares$naa))) %>% 
#   mutate(Weight=as.numeric(unlist(vpares$input$dat$waa)),
#          Maturity=as.numeric(unlist(vpares$input$dat$maa)))
# 

obs_w = as.matrix(waa_dat)[,1:3]
naa = as.matrix(vpares$naa)
minAge = 0
maxAgePlusGroup = 1

logwaa = as.matrix(log(vpares$input$dat$waa))

w0_dat = waa_dat %>% filter(Age==min(Age))
beta_w0 = log(mean(w0_dat$Weight))

A = waa_dat$Age %>% max

# waa_prev = sapply(1:nrow(waa_dat), function(i) {
#   if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
#     value <- NA
#   } else {
#     if (waa_dat$Age[i]<A) {
#       value <- as.numeric(
#         vpares$input$dat$waa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
#     } else {
#       value1 <- as.numeric(
#         vpares$input$dat$waa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
#       value2 <- as.numeric(
#         vpares$input$dat$waa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)]) # plus group
#       n1 <- as.numeric(
#         vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
#       n2 <- as.numeric(
#         vpares$naa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)]) # plus group
#       value <- (value1*n1+value2*n2)/(n1+n2)
#     }
#   }
#   value
# })
# 
# waa_dat = waa_dat %>%
#   mutate(Weight_prev=waa_prev) %>%
#   mutate(logN_prev=log(num_prev))
# %>%
#   # mutate(YearClass=Year-Age) %>%
#   # mutate(interact_norm = Weight_prev*Number_prev,
#   #        interact_log = Weight_prev*logN_prev) %>%
#   filter(Year<max(Year))

waa_dat2 = na.omit(waa_dat)

modw = glm(Weight~1+Weight_prev,data=waa_dat2,family=Gamma("identity"))
alpha_w = as.numeric(coef(modw)[1])
rho_w = as.numeric(coef(modw)[2])

logCV_w <- log(sqrt(summary(modw)$dispersion))



### growth dynamics with a state-space model 

maa = as.matrix(vpares$input$dat$maa)

maa_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    value <- as.numeric(
      vpares$input$dat$maa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
  }
  value
})

# picking up the same cohort(同じコホートにおける前年の資源尾数)
cohort_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if(waa_dat$Age[i]==0) {
      value <- NA
    } else {
      # if (waa_dat$Age[i]<A) {
      value <- as.numeric(vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      # } else {
      # }
    }
  }
  value
})

# 前年の同じコホートと1歳上のコホートの資源尾数の合計
cohort_plus = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if(waa_dat$Age[i]==0) {
      value <- NA
    } else {
      # if (waa_dat$Age[i]<A) {
      value <- as.numeric(vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      value <- value + as.numeric(vpares$naa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)])
      # } else {
      # }
    }
  }
  value
})


maa_dat = waa_dat %>% mutate(Maturity_prev = maa_prev,
                             Cohort_prev=cohort_prev,
                             Cohort_plus=cohort_prev)

delta_mat = function(Maturity,Maturity_prev) (Maturity-Maturity_prev)/(1-Maturity_prev)

maa_dat = maa_dat %>% 
  filter(Age>0 & Age<4) %>%
  # filter(Year > min(Year)) %>%
  mutate(y = delta_mat(Maturity,Maturity_prev))

obs_g = as.matrix(maa_dat %>% dplyr::select(Age,Year,Maturity))

maa_dat2 = maa_dat %>% na.omit() %>% mutate(Age = factor(Age))

mod_mat = glmmTMB(y~Age + Cohort_plus + Age:Cohort_plus + 1,data=maa_dat2,family=ordbeta)
mod_mat0 = glmmTMB(y~Age +  0,data=maa_dat2,family=ordbeta)
summary(mod_mat0)

mod_mat0$sdr

alpha_g = mod_mat0$fit$par[names(mod_mat0$fit$par)=="beta"] %>% as.numeric()
psi = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="psi"])
logdisp = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="betad"])[1]

ssb = as.numeric(colSums(vpares$ssb))

compile("growss.cpp")
dyn.load(dynlib("growss"))
# dyn.unload(dynlib("growss"))


g_fix = c(0,-1,-2,-3,1,1,1)
tmbdata0 = list(obs_w=obs_w,maa=maa,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
                dist_wobs=1,scale_number=1000,scale=1000,g_fix = c(0,-1,-2,-3,1,1,1))

# tmbdata0 = list(obs_w=obs_w,maa=maa,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
#                 dist_wobs=0,scale_number=1000,scale=1000,g_fix = c(0,-1,-2,-3,1,1,1),
#                 ssb=ssb)

# pars0=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#           omicron=log(0.1),logCV_w=logCV_w,alpha_g=alpha_g,psi=psi,logdisp=logdisp)

pars0=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
           omicron=log(0.1),logCV_w=logCV_w,alpha_g=alpha_g,psi=psi,logdisp=logdisp,beta_g=rep(0,3))

# pars0=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#            omicron=log(0.1),logCV_w=logCV_w)

obj0 = TMB::MakeADFun(tmbdata0,pars0,random="logwaa",DLL="growss")
# obj0 = TMB::MakeADFun(tmbdata0,pars0,random="logwaa",DLL="growss",
#                       map=list(beta_g=rep(factor(NA),3)))

opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
opt0

rep0 = sdreport(obj0,bias.correct = TRUE)
rep0

opt0$objective
# lognormal
# > opt0$objective
# [1] 2039.139

# gamma
# opt0$objective
# [1] 2038.827

# gammaの方が負の対数尤度が少しだけ小さい

opt0$par[names(opt0$par)=="psi"]
opt0$par[names(opt0$par)=="logdisp"]

REP0 = obj0$env$report()
REP0$maa_pred
matplot(cbind(ssb/1000,REP0$ssb))



# state space model for weight dynamics
compile("gloss_orig.cpp")
dyn.load(dynlib("gloss_orig"))

ssb = as.numeric(colSums(vpares$ssb))

# lognormal distribution
tmbdata1 = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
               dist_wobs=0,scale_number=1000,scale=1000,ssb=ssb)

pars=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
          iota=log(0.1),logCV_w=logCV_w)

obj1 = TMB::MakeADFun(tmbdata1,pars,random="logwaa",DLL="gloss_orig")

opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
opt1

rep1 = sdreport(obj1,bias.correct = TRUE)
rep1

# gamma distribution
tmbdata2 = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
                dist_wobs=1,scale_number=1000,scale=1000,ssb=ssb)

# pars=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#           iota=log(0.1),logCV_w=logCV_w)

obj2 = TMB::MakeADFun(tmbdata2,pars,random="logwaa",DLL="gloss_orig")

opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
opt2

rep2 = sdreport(obj2,bias.correct = TRUE)

rep1
rep2

c(opt1$objective,opt2$objective) #gamma distributionの方が負の対数尤度若干低い


# 密度効果なし

pars0=list(logwaa=logwaa,beta_w0=c(beta_w0),alpha_w=c(alpha_w),rho_w=rho_w,
          iota=log(0.1),logCV_w=logCV_w)


# pars=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#           iota=log(0.1),logCV_w=logCV_w)

obj0 = TMB::MakeADFun(tmbdata2,pars0,random="logwaa",DLL="gloss_orig")

opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
opt0

rep0 = sdreport(obj0,bias.correct = TRUE)
rep0
# 観測誤差がほぼ0になる

opt0

c(opt1$objective,opt2$objective,opt0$objective)

# AIC
2*c(opt1$objective,opt2$objective,opt0$objective)+
  2*c(length(opt1$par),length(opt2$par),length(opt0$par))


parList = obj2$env$parList()
names(parList)

parList[["beta_w0"]]
coef(mod_w0)
c(parList[["alpha_w"]],parList[["rho_w"]])

coef(mod_w_growth)
# 
# 
exp(parList[["logCV_w"]])
# 
exp(parList[["logwaa"]])
obj$report()[["sd_w"]]
obj$report()[["shape"]]
# どちらも0.1程度

# matrix(rep$unbiased$value,nrow=7)
# waa_obs

# w0（0歳魚の体重）をSSBでモデリングすると観測誤差も推定できる

dyn.unload(dynlib("gloss_orig"))

## try maturity modeling ----


maa_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    value <- as.numeric(
      vpares$input$dat$maa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
  }
  value
})

# picking up the same cohort(同じコホートにおける前年の資源尾数)
cohort_prev = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if(waa_dat$Age[i]==0) {
      value <- NA
    } else {
      # if (waa_dat$Age[i]<A) {
      value <- as.numeric(vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      # } else {
      # }
    }
  }
  value
})

# 前年の同じコホートと1歳上のコホートの資源尾数の合計
cohort_plus = sapply(1:nrow(waa_dat), function(i) {
  if (waa_dat$Year[i]==min(waa_dat$Year)) {
    value <- NA
  } else {
    if(waa_dat$Age[i]==0) {
      value <- NA
    } else {
      # if (waa_dat$Age[i]<A) {
      value <- as.numeric(vpares$naa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      value <- value + as.numeric(vpares$naa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)])
      # } else {
      # }
    }
  }
  value
})


maa_dat = waa_dat %>% mutate(Maturity_prev = maa_prev,
                             Cohort_prev=cohort_prev,
                             Cohort_plus=cohort_prev)

delta_mat = function(Maturity,Maturity_prev) (Maturity-Maturity_prev)/(1-Maturity_prev)

maa_dat = maa_dat %>% 
  filter(Age>0 & Age<4) %>%
  # filter(Year > min(Year)) %>%
  mutate(y = delta_mat(Maturity,Maturity_prev))

obs_g = as.matrix(maa_dat %>% dplyr::select(Age,Year,Maturity))

maa_dat2 = maa_dat %>% na.omit() %>% mutate(Age = factor(Age))

mod_mat = glmmTMB(y~Age + Cohort_plus + Age:Cohort_plus + 1,data=maa_dat2,family=ordbeta)
mod_mat0 = glmmTMB(y~Age +  0,data=maa_dat2,family=ordbeta)


temp = glmmTMB(y~1,data=maa_dat2,family=ordbeta)
temp$sdr

predict(temp,type="response") %>% mean
maa_dat2$y %>% mean

temp2 = glmmTMB(y~1,data=filter(maa_dat2,y>0 & y<1),family=beta_family)

predict(temp2,type="response") %>% mean
filter(maa_dat2,y>0 & y<1)$y %>% mean

mod_mat1 = glmmTMB(y~Age +  0,data=maa_dat2,family=ordbeta("probit"))

model_mat = glmmTMB(y ~ Age + Cohort_plus + Age:Cohort_plus + 0,
                    data=maa_dat2,family=ordbeta)

model_mat

# mod_mat1 = glmmTMB(y~Age +  0,data=maa_dat,family=ordbeta,dispformula = ~Age)
# mod_mat1$sdr

summary(mod_mat0)
summary(mod_mat1)

mod_mat0$obj$env$data
mod_mat1$obj$env$data

mod_mat0$obj$env$data$Xd

mean(maa_dat$y==0,na.rm=T)

summary(predict(mod_mat,type="response"))


logit = function(x) log(x/(1-x))
# logit = Vectorize(logit)

logit(c(0.1,0.2))

invlogit = function(x) 1/(1+exp(-x))
# invlogit = Vectorize(invlogit)
# plot(invlogit(-20:20))


compile("ordbeta.cpp")
dyn.load(dynlib("ordbeta"))


tmbdata = list(obs_g=obs_g)

alpha_g = mod_mat0$fit$par[names(mod_mat0$fit$par)=="beta"] %>% as.numeric()
psi = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="psi"])
logdisp = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="betad"])[1]

maa = vpares$input$dat$maa

logitmaa = matrix(
  rep(logit(rowMeans(maa[2:4,])),ncol(logwaa)),nrow=3)

naa0 = vpares$naa/1000

naa = naa0[-1,]+naa0[-nrow(naa0),] %>% as.matrix()

naa=as.matrix(naa[1:3,])
# cohort_plus

class(obs_g)
class(naa)

tmbdata = list(obs_g=obs_g,naa=naa)

# as.matrix(tmbdata$naa)

p=0.3
disp=10000
hist(rbeta(1000,p*disp,(1-p)*disp))

pars = list(alpha=alpha_g,logdisp=logdisp,psi=psi,omicron=log(1000),
            logitmaa=logitmaa,beta=rep(0,length(alpha_g)))

pars = list(alpha=alpha_g,logdisp=log(1000),psi=psi,omicron=1,
            logitmaa=logitmaa,beta=rep(0,length(alpha_g)))


map=list()
# map = list(omicron=factor(NA))
map = list(logdisp=factor(NA))

# map = list(omicron=factor(NA),psi=rep(factor(NA),2))

obj = TMB::MakeADFun(tmbdata,pars,random=c("logitmaa"),DLL="ordbeta",
                     map=map)

# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss_orig")

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

rep = sdreport(obj,bias.correct = TRUE)
rep


alpha_g
alpha_g = c(-3.5,-0.5,1)
omicron=1
logdisp=1
disp=exp(logdisp)
# mu = c(0.1)
# s1 = mu*disp
# s2 = (1-mu)*disp
# r=rnorm(n,0,1)
# v=pnorm(r,0,1) 
# w=qbeta(v,s1,s2)
maa_obs <- maa_true  <- matrix(0,nrow=3,ncol=50)
maa_true[,1] <- c(0.15,0.75,0.95)
set.seed(1)
for(j in 2:ncol(maa_true)){
  for(i in 1:3) {
    if(i==1){
      tmp = invlogit(alpha_g[i])
    } else{
      tmp = maa_true[i-1,j-1]+invlogit(alpha_g[i])*(1-maa_true[i-1,j-1])
    }
    s1 = tmp*exp(omicron)
    s2 = (1-tmp)*exp(omicron)
    r = rnorm(1,0,1)
    v = pnorm(r,0,1)
    w = qbeta(v,s1,s2)
    maa_true[i,j] <- w
  }
}
set.seed(2)
for(j in 1:ncol(maa_true)){
  for(i in 1:3) {
    maa_obs[i,j] <- rbeta(1,maa_true[i,j]*disp,(1-maa_true[i,j])*disp)
    maa_obs[i,j] <- round(maa_obs[i,j],2)
  }
}

maa_obs

# compile("ordbeta.cpp")
# dyn.load(dynlib("ordbeta"))

obs_g  = data.frame(age=rep(1:nrow(maa_true)-1,ncol(maa_true)),
                    year=rep(1:ncol(maa_true)-1,each=nrow(maa_true)),
                    maturity=as.numeric(maa_obs)) %>%
  # mutate(maturity = y_trans(maturity,n())) %>%
  as.matrix()

logitmaa = log(maa_true)/(1-maa_true)
tmbdata=list(obs_g=obs_g)
pars=list(alpha=alpha_g,logdisp=logdisp,psi=c(0,0),omicron=omicron,logitmaa=logitmaa)

# map=list(psi=rep(factor(NA),2))
# map=list(alpha=rep(factor(NA),length(alpha_g)),logdisp=factor(NA),omicron=factor(NA),
#          psi=rep(factor(NA),2))
# map=list(logdisp=factor(NA),omicron=factor(NA),
#          psi=rep(factor(NA),2))

map=list()
# map=list(omicron=factor(NA))
obj = TMB::MakeADFun(tmbdata,pars,random=c("logitmaa"),DLL="ordbeta",
                     map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

opt
rep = sdreport(obj)
rep

dyn.unload(dynlib("ordbeta"))







compile("gloss.cpp")
dyn.load(dynlib("gloss"))

# dyn.unload(dynlib("gloss"))
# vpares$input$dat$maa

compile("gloss_orig.cpp")
dyn.load(dynlib("gloss_orig"))


g_fix = c(0,-1,-2,-3,1,1,1)

# g_fix = c(0,-1,-1,-1,1,1,1)

# dim(obs_g)

maa = vpares$input$dat$maa %>% as.matrix()

# tmbdata = list(obs_w=obs_w,maa=maa,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
#                dist_wobs=0,scale_number=1000,g_fix=g_fix,scale=1000,ssb=ssb)

tmbdata = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
               dist_wobs=1,scale_number=1000,g_fix=g_fix,scale=1000,ssb=ssb)



# coef(summary(mod_mat0))
alpha_g = mod_mat0$fit$par[names(mod_mat0$fit$par)=="beta"] %>% as.numeric()
psi = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="psi"])
omicron = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="betad"])[1]

# logit_deltag=matrix(rep(alpha_g,ncol(logwaa)),nrow=3)

maa = vpares$input$dat$maa
logitmaa = matrix(
  rep(logit(rowMeans(maa[2:4,])),ncol(logwaa)),nrow=3)

# logitmaa = matrix(
#   0,ncol=ncol(logwaa),nrow=3)
# 

# logit_deltag=matrix(predict(mod_mat0),nrow=3)


# dim(logwaa)
# dim(logitg)
# dim(naa)

pars=list(logwaa=logwaa,
          beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
          iota=log(0.1),omicron=omicron,logCV_w=logCV_w,alpha_g=alpha_g,
          psi=psi)


pars=list(logwaa=logwaa,
          beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
          iota=log(0.1),logCV_w=logCV_w)

pars=list(logwaa=logwaa,
          beta_w0=c(beta_w0,0)[1],alpha_w=c(alpha_w,0)[1],rho_w=rho_w,
          iota=log(0.1),logCV_w=logCV_w)


# 
# pars=list(logwaa=logwaa,logitmaa=logitmaa,
#           beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#           iota=log(0.1),omicron=log(5),logCV_w=logCV_w,alpha_g=alpha_g,
#           psi=psi,logdisp=logdisp-1)

# pars=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#           iota=log(0.1),logCV_w=logCV_w,alpha_g=alpha_g,
#           psi=psi,logdisp=logdisp)
# 

# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa","logitmaa"),DLL="gloss")
obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss")

obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss_orig")

# map = list()

# pars$maa
# 
# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss",
#                      map=map)
# 
# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa"),DLL="gloss")

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

rep = sdreport(obj,bias.correct = TRUE)
rep

obj$env$report()

psi



dyn.unload(dynlib("gloss"))


# tray simple ordbeta model

compile("ordbeta.cpp")
dyn.load(dynlib("ordbeta"))

tmp = maa_dat %>%
  mutate(x=Age-min(Age)) %>%
  # filter(0<y & y<1) %>% 
  dplyr::select(y,x) %>%
  na.omit()

y_trans = function(y0,N) (y0*(N-1)+0.5)/N

y = y_trans(tmp$y,nrow(tmp))
x = rep(0,length(y))
x = tmp$x

tmp
# res = glmmTMB(y~1,data=tmp,family=beta_family)
res = glmmTMB(y~0+factor(x),data=tmp,family=ordbeta)
res
res$sdr


alpha = res$fit$par[names(res$fit$par)=="beta"]
logdisp = res$fit$par[names(res$fit$par)=="betad"]
# psi = res$fit$par[names(res$fit$par)=="psi"]

tmbdata=list(y=y,x=x)
# pars=list(alpha=alpha,logdisp=logdisp,psi=psi)
# pars=list(alpha=alpha,logdisp=logdisp)

obj = TMB::MakeADFun(tmbdata,pars,DLL="ordbeta")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt
rep = sdreport(obj)
rep

res$sdr

dyn.unload(dynlib("ordbeta"))


n=100

logdisp=2
disp=exp(logdisp)
# disp=100
# sd=0.1
mu = c(0.1)
s1 = mu*disp
s2 = (1-mu)*disp
y = rbeta(n,s1,s2)
hist(y)

logdisp=2
disp=exp(logdisp)
mu = c(0.1)
s1 = mu*disp
s2 = (1-mu)*disp
r=rnorm(n,0,1)
v=pnorm(r,0,1) 
w=qbeta(v,s1,s2)

hist(w)


v=pbeta(w,s1,s2)
r=qnorm(v,0,1)

r=rnorm(n,0,sd)
y=logit(mu)+r
hist(invlogit(y))



# mu=pnorm(r,0,1)
# logdisp=2
# disp=exp(logdisp)
# disp=1000
s1=mu*disp
s2=(1-mu)*disp
w = qbeta(mu,s1,s2)
hist(w)

i<-1

alpha_g
alpha_g = c(-3.5,-0.5,1)
omicron=1
logdisp=1
disp=exp(logdisp)
# mu = c(0.1)
# s1 = mu*disp
# s2 = (1-mu)*disp
# r=rnorm(n,0,1)
# v=pnorm(r,0,1) 
# w=qbeta(v,s1,s2)
maa_obs <- maa_true  <- matrix(0,nrow=3,ncol=50)
maa_true[,1] <- c(0.15,0.75,0.95)
set.seed(1)
for(j in 2:ncol(maa_true)){
  for(i in 1:3) {
    if(i==1){
      tmp = invlogit(alpha_g[i])
    } else{
      tmp = maa_true[i-1,j-1]+invlogit(alpha_g[i])*(1-maa_true[i-1,j-1])
    }
    s1 = tmp*exp(omicron)
    s2 = (1-tmp)*exp(omicron)
    r = rnorm(1,0,1)
    v = pnorm(r,0,1)
    w = qbeta(v,s1,s2)
    maa_true[i,j] <- w
  }
}
set.seed(2)
for(j in 1:ncol(maa_true)){
  for(i in 1:3) {
    maa_obs[i,j] <- rbeta(1,maa_true[i,j]*disp,(1-maa_true[i,j])*disp)
    maa_obs[i,j] <- round(maa_obs[i,j],2)
  }
}

maa_obs

compile("ordbeta.cpp")
dyn.load(dynlib("ordbeta"))

obs_g  = data.frame(age=rep(1:nrow(maa_true)-1,ncol(maa_true)),
                    year=rep(1:ncol(maa_true)-1,each=nrow(maa_true)),
                    maturity=as.numeric(maa_obs)) %>%
  mutate(maturity = y_trans(maturity,n())) %>%
  as.matrix()

logitmaa = log(maa_true)/(1-maa_true)
tmbdata=list(obs_g=obs_g)
# pars=list(alpha=alpha_g,logdisp=logdisp,psi=c(0,0),omicron=omicron,logitmaa=logitmaa)

logsd = log(sd(as.numeric(maa_obs-maa_true)))

# pars=list(alpha=alpha_g,logdisp=logdisp,omicron=omicron,logitmaa=logitmaa)
# pars=list(alpha=alpha_g,logdisp=logdisp,psi=c(0,0),omicron=log(0.05),logitmaa=logitmaa)
pars=list(alpha=alpha_g,logdisp=logdisp,omicron=log(0.05),logitmaa=logitmaa)

# map=list(psi=rep(factor(NA),2))
# map=list(alpha=rep(factor(NA),length(alpha_g)),logdisp=factor(NA),omicron=factor(NA),
#          psi=rep(factor(NA),2))
# map=list(logdisp=factor(NA),omicron=factor(NA),
#          psi=rep(factor(NA),2))

map=list()
# map=list(omicron=factor(NA))
obj = TMB::MakeADFun(tmbdata,pars,random=c("logitmaa"),DLL="ordbeta",
                     map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

opt
rep = sdreport(obj)
rep

dyn.unload(dynlib("ordbeta"))


