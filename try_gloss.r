

### set workind directory ---- 

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

# compile("gloss.cpp")
# dyn.load(dynlib("gloss"))

ssb = as.numeric(colSums(vpares$ssb))

# tmbdata = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
#                dist_wobs=0,scale_number=1000,scale=1000,ssb=ssb)
# 
# tmbdata = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
#                dist_wobs=1,scale_number=1000,scale=1000,ssb=ssb)


# tmbdata = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
#                dist_wobs=1)
# 
# pars=list(logwaa=logwaa,beta_w0=beta_w0,alpha_w=alpha_w,rho_w=rho_w,
#           iota=log(0.1),logCV_w=logCV_w)

# pars=list(logwaa=logwaa,beta_w0=beta_w0,alpha_w=c(alpha_w,0),rho_w=rho_w,
#           iota=log(0.1),logCV_w=logCV_w)
# 
# pars=list(logwaa=logwaa,beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
#           iota=log(0.1),logCV_w=logCV_w)


# tmbdata$obs_w

# obj = TMB::MakeADFun(tmbdata,pars,random="logwaa",DLL="gloss")
# 
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# opt
# 
# rep = sdreport(obj,bias.correct = TRUE)
# rep
# 
# 
# parList = obj$env$parList()
# names(parList)
# 
# parList[["beta_w0"]]
# coef(mod_w0)
# 
# c(parList[["alpha_w"]],parList[["rho_w"]])
# coef(mod_w_growth)
# 
# 
# exp(parList[["logCV_w"]])
# 
# exp(parList[["logwaa"]])
# obj$report()[["sd_w"]]
# obj$report()[["shape"]]

# どちらも0.1程度

# matrix(rep$unbiased$value,nrow=7)
# waa_obs

# 観測誤差がほぼ0になる
# w0（0歳魚の体重）をSSBでモデリングすると観測誤差も推定できる



# try maturity modeling

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

mod_mat1 = glmmTMB(y~Age +  0,data=maa_dat2,family=ordbeta("probit"))

# mod_mat1 = glmmTMB(y~Age +  0,data=maa_dat,family=ordbeta,dispformula = ~Age)
# mod_mat1$sdr

summary(mod_mat0)
summary(mod_mat1)

mod_mat0$obj$env$data
mod_mat1$obj$env$data

mod_mat0$obj$env$data$Xd

mean(maa_dat$y==0)

summary(predict(mod_mat,type="response"))

# invlogit(-1.9925119)

# mod_mat0$sdr
# 
# parl = mod_mat0$obj$env$parList()
# names(parl)
# parl[["psi"]]
# parl[["betad"]] 
# 
# parl[["betad"]] %>% exp
# # beta分布のdispersion parameter
# 
# summary(mod_mat0)
# 
# mod_mat0$obj$env$map
# 
# invlogit = function(x) 1/(1+exp(-x))
# invlogit = Vectorize(invlogit)
# plot(invlogit(-20:20))

compile("gloss.cpp")
dyn.load(dynlib("gloss"))

# vpares$input$dat$maa

g_fix = c(0,-1,-2,-3,1,1,1)

# dim(obs_g)

tmbdata = list(obs_w=obs_w,obs_g=obs_g,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
               dist_wobs=0,scale_number=1000,g_fix=g_fix,scale=1000,ssb=ssb)


# coef(summary(mod_mat0))
alpha_g = mod_mat0$fit$par[names(mod_mat0$fit$par)=="beta"] %>% as.numeric()
psi = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="psi"])
logdisp = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="betad"])[1]

logit_deltag=matrix(rep(alpha_g,ncol(logwaa)),nrow=3)

# logit_deltag=matrix(predict(mod_mat0),nrow=3)



dim(logwaa)
dim(logit_deltag)
dim(naa)

pars=list(logwaa=logwaa,logit_deltag=logit_deltag,
          beta_w0=c(beta_w0,0),alpha_w=c(alpha_w,0),rho_w=rho_w,
          iota=log(0.1),omicron=log(0.1),logCV_w=logCV_w,alpha_g=alpha_g,
          psi=psi,logdisp=logdisp)

obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa","logit_deltag"),DLL="gloss")

# map = list(omicron=factor(NA),alpha_g=rep(factor(NA),length(alpha_g)),
#            psi=rep(factor(NA),length(psi)),logdisp=factor(NA))
# 
# obj = TMB::MakeADFun(tmbdata,pars,random=c("logwaa","logit_deltag"),DLL="gloss",
#                      map=map)

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

rep = sdreport(obj,bias.correct = TRUE)
rep

