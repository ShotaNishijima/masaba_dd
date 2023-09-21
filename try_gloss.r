

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

scale_num = A = waa_dat$Age %>% max

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


waa_dat = expand.grid(Age=as.numeric(rownames(vpares$naa)),
                      Year=as.numeric(colnames(vpares$naa))) %>% 
  mutate(Weight=as.numeric(unlist(vpares$input$dat$waa)),
         Maturity=as.numeric(unlist(vpares$input$dat$maa)))




obs_w = as.matrix(waa_dat)[,1:3]
naa = as.matrix(vpares$naa)
minAge = 0
maxAgePlusGroup = 1

logwaa = as.matrix(log(vpares$input$dat$waa))

w0_dat = waa_dat %>% filter(Age==min(Age))
beta_w0 = log(mean(w0_dat$Weight))

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

waa_dat = waa_dat %>%
  mutate(Weight_prev=waa_prev) %>%
  mutate(logN_prev=log(num_prev))
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

compile("gloss.cpp")
dyn.load(dynlib("gloss"))

tmbdata = list(obs_w=obs_w,naa=naa,minAge=minAge,maxAgePlusGroup=maxAgePlusGroup,
               dist_wobs=0) # lognormal distribusion

pars=list(logwaa=logwaa,beta_w0=beta_w0,alpha_w=alpha_w,rho_w=rho_w,
          iota=log(0.1),logCV_w=logCV_w)

tmbdata$obs_w


obj = TMB::MakeADFun(tmbdata,pars,random="logwaa",DLL="gloss")

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

rep = sdreport(obj)
rep

# 観測誤差がほぼ0になる




rep$value

parList = obj$env$parList()
names(parList)

exp(parList[["logwaa"]])
waa
obj$report()[["sd_w"]]

waa_obs
