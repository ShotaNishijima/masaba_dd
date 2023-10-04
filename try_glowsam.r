
### set working directory ---- 

gitdir = "~/git/masaba_dd"
setwd(gitdir)

savedir = "res"

if (!file.exists(savedir)) dir.create(savedir)

savename = function(x) paste0(savedir,"/",x)


library(frasyr)
# library(frasam)
devtools::load_all("~/git/frasam")
library(magrittr)
library(tidyverse)

vpares = get(load("data/vpa_masaba_P2022.rda")) # assessment result in 2022

use_sam_tmb(overwrite=TRUE)

dat = vpares$input$dat

dat$index

#最新年の加入indexが抜けているので元データをとってくる
cpue = read.csv("~/git/masaba2022/Data/cpue_2022_final.csv",row.names=1)
dat$index[,ncol(dat$index)] <- cpue[,ncol(cpue)]


# かなりシンプルなモデル
samres0 <- sam(
  dat, #rvpaと同じデータ
  # dat2,
  last.catch.zero = TRUE,
  abund = c("N","N","SSB","SSB","N"),
  sel.def="max",
  min.age=c(0,0,0,0,0,1),
  max.age = c(0,0,6,6,6,1),
  rec.age = 0,
  b.est=TRUE,
  b.fix=c(NA,NA,1,1,NA), #VPAの設定と同じ
  SR = "BH", 
  varC = c(0,1,2,2,3,0,0),
  varF = c(0,0,1,1,1,1,1),
  varN = c(0,1,1,1,1,1,1),
  varN.fix=c(NA,1e-4), #New!
  rho.mode=3,
  bias.correct = TRUE, #IMPORTANT
  get.random.vcov = FALSE,
  # remove.Fprocess.year = 2011, # IMPORTANT
  scale=1000,
  gamma=100000,
  use.index=1:5,
  est.method = "ls",
  add_random = c("rec_logb","rec_loga")[1]
  # index.key=c(1,2,3,4,5,5),
  # index.key=c(0,0,1,2,0) #indexのsigmaの設定
  # index.b.key = c(0,1,2,3,4,4) # 
)

samres0$rep$par.random %>% names %>% unique

samres0$rep$par.random[names(samres0$rep$par.random)=="rec_logb"]

# samres0$rec.par

# samres0$rep

matplot(cbind(colSums(samres3$ssb),colSums(samres0$ssb)))

samres0$b
samres0$sigma

samres1 <- sam(
  dat, #rvpaと同じデータ
  # dat2,
  last.catch.zero = TRUE,
  abund = c("N","N","SSB","SSB","N"),
  sel.def="max",
  min.age=c(0,0,0,0,0,1),
  max.age = c(0,0,6,6,6,1),
  rec.age = 0,
  b.est=TRUE,
  b.fix=c(1.8,1.8,1,1,1.8), #VPAの設定と同じ
  SR = "BH", 
  varC = c(0,1,2,2,3,4,4),
  varF = c(0,0,1,1,1,1,1),
  varN = c(0,1,1,1,1,1,1),
  varN.fix=c(NA,1e-4), #New!
  rho.mode=3,
  bias.correct = TRUE, #IMPORTANT
  get.random.vcov = FALSE,
  # remove.Fprocess.year = 2011, # IMPORTANT
  scale=1000,
  gamma=100000,
  use.index=1:5,
  # index.key=c(1,2,3,4,5,5),
  index.key=c(0,0,1,2,0) #indexのsigmaの設定
  # index.b.key = c(0,1,2,3,4,4) # 
)

input = samres1$input
input$b.fix <- c(NA,NA,1,1,NA)
samres2 = do.call(sam,input)
samres2$b
samres2$sigma
samres2$aic

input = samres2$input
input$index.key <- 0:4
samres3 = do.call(sam,input)
samres3$b
samres3$sigma
samres2$sigma
c(samres1$aic,samres2$aic,samres3$aic)

input = samres2$input
input$index.b.key <- c(0,0,1,1,2)
input$b.fix <- c(NA,1,NA)
samres4 = do.call(sam,input)
samres4$sigma
samres4$b
c(samres1$aic,samres2$aic,samres3$aic,samres4$aic)


input = samres3$input
input$tmbdata = samres3$data
input$map = samres3$map
samres5 = do.call(sam,input)

# samres5$sigma.logN
# 
c(samres3$loglik,samres5$loglik)

dyn.unload(dynlib("sam"))


### exploring new sam ----


waa_dat = expand.grid(Age=as.numeric(rownames(dat$waa)),
                      Year=as.numeric(colnames(dat$waa))) %>% 
  mutate(Weight=as.numeric(unlist(dat$waa)),
         Maturity=as.numeric(unlist(dat$maa)))

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

samres3$naa

waa_dat = waa_dat %>%
  mutate(Weight_prev=waa_prev) %>%
  mutate(YearClass=Year-Age) 

modw = glm(Weight~1+Weight_prev,data=waa_dat2,family=Gamma("identity"))
alpha_w = as.numeric(coef(modw)[1])
rho_w = as.numeric(coef(modw)[2])
logCV_w <- log(sqrt(summary(modw)$dispersion))

beta_w0 = waa_dat %>% filter(Age==min(Age)) %>% pull(.,Weight) %>% mean %>% log

modw2 = glm(log(Weight)~1+Weight_prev,data=waa_dat2)
sd(modw2$residuals)

tmp = waa_dat %>% filter(Age==min(Age)) %>% pull(.,Weight) %>% log %>% sd

logCV_w <- c(log(tmp),logCV_w)

use_sam_tmb(TmbFile = "glowsam",overwrite=TRUE)

RES = samres3
RES = samres0
RES$b
tmbdata = RES$data
map0 = RES$map
obj0 = RES$obj
# obj0$env$.random

opt0 = RES$opt
rep0 = RES$rep
  
parList0 = RES$par_list
# tmbdata$model_weight_maturity = c(-1,-1)
tmbdata$model_weight_maturity = c(1,-1)
tmbdata$scale_number = 1000
tmbdata$dist_wobs = 0 #lognormal(0) or gamma(1)
weight_mat = as.matrix(dat$waa)
weight_mat[] <- 1
tmbdata$weight_weight <- weight_mat
tmbdata$maturity_weight <- weight_mat

parList <- parList0
parList$logwaa <- as.matrix(log(dat$waa))
parList$beta_w0 = beta_w0
# parList$beta_w0 = c(beta_w0,0)
parList$alpha_w = c(alpha_w)
# parList$alpha_w = c(alpha_w,0)
parList$rho_w = rho_w
parList$omicron = log(0.1)
parList$logCV_w <- logCV_w

# parList$logwaa %>% class
# class(tmbdata$weight_weight)

map <- map0

map$rec_logb <- factor(NA)

# map$logwaa <- factor(matrix(NA,ncol=ncol(dat$waa),nrow=nrow(dat$waa)))
# map$alpha_w <- rep(factor(NA),length(parList$alpha_w))
# map$beta_w0 <- rep(factor(NA),length(parList$beta_w0))
# map$rho_w <- rep(factor(NA),length(parList$rho_w))
# map$omicron <- factor(NA)
# map$logCV_w <- rep(factor(NA),2)


# obj = MakeADFun(data=tmbdata,parameters=parList,map=map,
#                 random=c("U"),DLL="glowsam")

obj = MakeADFun(data=tmbdata,parameters=parList,map=map,
                random=c("U","logwaa"),DLL="glowsam")

# obj = MakeADFun(data=tmbdata,parameters=parList,map=map,
#                 random=c("U","logwaa","rec_logb"),DLL="glowsam")

opt <- nlminb(obj$par, obj$fn, obj$gr,control=list(eval.max=1000,iter.max=1000))
opt
rep <- sdreport(obj,bias.correct = TRUE)

rep$value %>% names %>% unique

SSB= rep$unbiased$value[names(rep$value)=="ssb"]
NAA = matrix(rep$unbiased$value[names(rep$value)=="exp_logN"],nrow=7)
WAA = t(matrix(rep$unbiased$value[names(rep$value)=="stockMeanWeight_true"],nrow=length(dat$waa)))
colnames(WAA) <- colnames(dat$waa)

WAA
waa0

R = matrix(rep$unbiased$value[names(rep$value)=="exp_logN"],nrow=7)[1,]
plot(R~SSB,type="l")


naa0 = samres0$naa

cbind(naa0[1,],R)

maa = dat$maa

# WAA = exp(matrix(rep$par.random[names(rep$par.random)=="logwaa"],nrow=7))


waa0 = dat$waa
colSums(naa0*waa0*maa)

colSums(NAA*WAA*maa)

plot(colSums(NAA*WAA*maa)/1000)
plot(SSB)

cbind(
  samres0$rep$unbiased$value[names(samres0$rep$value)=="ssb"],colSums(NAA*WAA*maa)
) %>% matplot


cbind(
  samres0$rep$unbiased$value[names(samres0$rep$value)=="ssb"],SSB)

colSums(samres0$ssb)/1000

names(rep$value) %>% unique

rep$par.random
names(rep$par.random) %>% unique

rep$par.random[names(rep$par.random)=="rec_logb"]

exp(-1.73)
exp(-2.27)

RES$sigma
RES$sigma.logC

opt
rep


rep$gradient.fixed

rep$gradient.fixed[22]
opt$par[22]

rep

tmbdata$stockMeanWeight
obj$env$report()[["stockMeanWeight_true"]]

# 密度効果を入れないと観測誤差のCVがもろ小さくなる

rep$sd[names(rep$value)=="stockMeanWeight_true"]


cbind(opt0$par,opt$par)


tmbdata$catchMeanWeight
tmbdata$stockMeanWeight
tmbdata$landMeanWeight



samres4$rec.par



