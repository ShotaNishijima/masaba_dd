
###############################################################################
### Estimating Population-level density dependence and MSY reference points ###
### by incorporating individual-level (post-recruit) density dependence #######
############################### Shota Nishijima ###############################

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

### read and handle data ----

# VPA estimates from the Japanese stock assessment in 2021
# vpares = get(load("data/vpa_masaba_P2021.rda"))
# vpares = get(load("data/res_vpa_CMP.rda")) # assessment result in 2021

vpares = get(load("data/vpa_masaba_P2022.rda")) # assessment result in 2022

range(colSums(vpares$ssb))
range(colSums(vpares$ssb))/1000

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

### modeling individual-level density dependence ----

## weight growth modeling ----

# Brody-growth coefficient and the von-Bertalanffy weight model
# https://cdn.inst-fs-iad-prod.inscloudgate.net/6c9d7269-b2b3-4885-8eaf-ce77b8700f70/24%20Delay-difference%20model.pdf?token=eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCIsImtpZCI6ImNkbiJ9.eyJyZXNvdXJjZSI6Ii82YzlkNzI2OS1iMmIzLTQ4ODUtOGVhZi1jZTc3Yjg3MDBmNzAvMjQlMjBEZWxheS1kaWZmZXJlbmNlJTIwbW9kZWwucGRmIiwidGVuYW50IjoiY2FudmFzIiwidXNlcl9pZCI6bnVsbCwiaWF0IjoxNjk1MjQ5MzY4LCJleHAiOjE2OTUzMzU3Njh9.wPWvRm1Ct7-q5K7vd18HENBOeM_t-eAaMglUhdvRpSp1GTRa6b-Aw-Lx1LCLEIjqoF-Lv1Vz0xdl13OIbVKiig&download=1&content_type=application%2Fpdf
# https://openmse.com/features-assessment-models/1-dd/

# w_{a+1} = alpha + rho*w_{a}
# rho = (w_{a+1}-w_inf/w_{a}-w_inf)
# alpha = w_inf*(1-rho)

head(waa_dat)
colnames(waa_dat)

waa_dat2 = na.omit(waa_dat)

colnames(waa_dat2)

(gg1 = ggplot(waa_dat2,aes(x=Number_prev,group=factor(Age),fill=factor(Age)),alpha=0.2)+
  geom_histogram(position="identity"))

(gg2 = ggplot(waa_dat2,aes(x=log(Number_prev),fill=factor(Age)),alpha=0.2)+
    geom_histogram(position="identity"))

(gg3 = ggplot(waa_dat2,aes(x=log(Cohort_prev),fill=factor(Age)),alpha=0.8)+
    geom_histogram(position="stack"))

(gg4 = ggplot(waa_dat2,aes(x=Cohort_prev,fill=factor(Age)),alpha=0.8)+
    geom_histogram(position="identity"))

waa_dat2 = na.omit(waa_dat) %>% group_by(Age) %>%
  mutate(logN_mean = mean(log(Number_prev)),
         logn_mean = mean(log(Cohort_prev)),
         logn_plus_mean = mean(log(Cohort_plus)),
         logN_sd = sd(log(Number_prev)),
         logn_sd = sd(log(Cohort_prev)),
         logn_plus_sd = sd(log(Cohort_plus))) %>%
  # ungroup() %>%
  # mutate(DF_N = log(Number_prev),
  #        DF_n = log10(Cohort_prev/n_mean),
  #        DF_n_plus = log10(Cohort_plus/n_plus_mean))
  mutate(DF_N = scale(log(Number_prev))[,1],
         DF_n = scale(log(Cohort_prev))[,1],
         DF_n_plus = scale(log(Cohort_plus))[,1]
  ) %>% ungroup

# plot(waa_dat2$DF_N,waa_dat3$DF_N)

# waa_dat3 = na.omit(waa_dat) %>% group_by(Age) %>%
#   mutate(N_mean = mean(Number_prev),
#          n_mean = mean(Cohort_prev),
#          n_plus_mean = mean(Cohort_plus)) %>%
#   ungroup() %>%
#   # mutate(DF_N = log(Number_prev),
#   #        DF_n = log10(Cohort_prev/n_mean),
#   #        DF_n_plus = log10(Cohort_plus/n_plus_mean))
#   mutate(DF_N = log10(Number_prev/N_mean),
#          DF_n = log10(Cohort_prev/n_mean),
#          DF_n_plus = log10(Cohort_plus/n_plus_mean))
# 

# 個体数をそのまま使うように変更(2022/04/09)
full_w_growth = glm(Weight~Weight_prev*Number_prev+Weight_prev*Cohort_prev+Weight_prev*Cohort_plus,
                    data=waa_dat2,family=Gamma("identity"))

# full_w_growth = glm(Weight~Weight_prev*DF_N+Weight_prev*DF_n+Weight_prev*DF_n_plus,
#                     data=waa_dat2,family=Gamma("identity"))

summary(full_w_growth)

# (tmp_aicc = AICc(mod_w_growth)) #3461.547

# getAllTerms(full_w_growth)

plot(waa_dat2$Cohort_plus~waa_dat2$Cohort_prev,log="xy")

dredge_w_growth = dredge(full_w_growth,fixed="Weight_prev",
                         subset=!("Number_prev" && "Cohort_prev") & 
                           !("Cohort_plus" && "Cohort_prev") & !("Number_prev" && "Cohort_plus"))

dredge_w_growth

head(dredge_w_growth,100)

# head(dredge_w_growth2,100)

save(dredge_w_growth,file=savename("dredge_w_growth.rda"))
# model.sel(dredge_w_growth)

model_sel_w_growth = model.sel(dredge_w_growth)

# model_sel_w_growth = model.sel(dredge_w_growth,beta="sd")

write.csv(as.data.frame(model_sel_w_growth),file=savename("AICc_table_w_growth.csv"))

mod_w_growth = get.models(dredge_w_growth,subset=1)[[1]]
mod_w_growth$formula
# Weight ~ Number_prev + 1 + Weight_prev
# alphaに前年の個体数の総数が影響する
summary(mod_w_growth)

AICc(mod_w_growth)

save(mod_w_growth,file=savename("model_w_growth.rda"))

## figure weight growth ----

base_size=12
point_size=1.5
path_size=1.2

waa_dat2$Number_prev %>% summary()
newdata_wg = expand.grid(Weight_prev=seq(min(waa_dat2$Weight_prev)*1,max(waa_dat2$Weight_prev)*1,length=100),
                         Number_prev=c(quantile(waa_dat2$Number_prev,probs=c(0.1,0.9)),mean(waa_dat2$Number_prev))) %>%
  as.data.frame()

# newdata_wg = expand.grid(Weight_prev=seq(min(waa_dat2$Weight_prev)*1,max(waa_dat2$Weight_prev)*1,length=100),
#                          DF_n_plus=c(quantile(waa_dat2$DF_n_plus,probs=c(0.1,0.9)),mean(waa_dat2$DF_n_plus))) %>%
#   as.data.frame()


# newdata_wg = expand.grid(Weight_prev=seq(min(waa_dat2$Weight_prev)*1,max(waa_dat2$Weight_prev)*1,length=100),
#                          DF_N=c(quantile(waa_dat2$DF_N,probs=c(0.1,0.9)),mean(waa_dat2$DF_N))) %>%
#   as.data.frame()


newdata_wg = newdata_wg %>% mutate(Weight=predict(mod_w_growth,newdata=newdata_wg))

newdata_wg$Number_prev %>% unique()


(g_wg = ggplot(data=NULL,aes(x=Weight_prev,y=Weight))+
  geom_point(data=waa_dat2,aes(colour=Number_prev),size=point_size)+
    # xlim(0,NA)+ylim(0,NA)+
    geom_path(data=newdata_wg,aes(colour=Number_prev,group=Number_prev),linewidth=path_size)+
   scale_colour_gradient(low="deepskyblue",high="sienna1",name=bquote(italic(N)))+
   theme_bw(base_size=base_size)+
        xlab(bquote(italic(w[t])))+ylab(bquote(italic(w[t+1])))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0.02, 0.02))
)


# (g_wg = ggplot(data=NULL,aes(x=Weight_prev,y=Weight))+
#     geom_point(data=waa_dat2,aes(colour=DF_n_plus),size=point_size)+
#     # xlim(0,NA)+ylim(0,NA)+
#     geom_path(data=newdata_wg,aes(colour=DF_n_plus,group=DF_n_plus),size=path_size)+
#     scale_colour_gradient(low="deepskyblue",high="sienna1",name=bquote('DF('*italic(n)*"+)"))+
#     theme_bw(base_size=base_size)+
#     xlab(bquote(italic(w[t])))+ylab(bquote(italic(w[t+1])))+
#     expand_limits(x = 0, y = 0)+
#     scale_x_continuous(expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0.02, 0.02))
#   # +geom_abline(data=NULL,aes(intercept=0,slope=1))
# )

# (g_wg = ggplot(data=NULL,aes(x=Weight_prev,y=Weight))+
#     geom_point(data=waa_dat2,aes(colour=DF_N),size=point_size)+
#     # xlim(0,NA)+ylim(0,NA)+
#     geom_path(data=newdata_wg,aes(colour=DF_N,group=DF_N),size=path_size)+
#     scale_colour_gradient(low="deepskyblue",high="sienna1",name=bquote('DF('*italic(N)*")"))+
#     theme_bw(base_size=base_size)+
#     xlab(bquote(italic(w[t])))+ylab(bquote(italic(w[t+1])))+
#     expand_limits(x = 0, y = 0)+
#     scale_x_continuous(expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0.02, 0.02))
#   # +geom_abline(data=NULL,aes(intercept=0,slope=1))
# )

ggsave(g_wg,filename=savename("weight_growth.png"),dpi=600,height=100,width=150,unit="mm")

save(g_wg,file=savename("weight_growth_graph.rda"))

## initial weight modeling ----

w0_dat = waa_dat %>% filter(Age==0 ) 

rec = as.numeric(vpares$naa[1,as.character(w0_dat$Year)])/1000
ssb = as.numeric(colSums(vpares$ssb[,as.character(w0_dat$Year-0)]))/1000
biom = as.numeric(colSums(vpares$baa[,as.character(w0_dat$Year-0)]))/1000
num = as.numeric(colSums(vpares$naa[,as.character(w0_dat$Year-0)]))/1000
ssn = as.numeric(colSums(vpares$naa[,as.character(w0_dat$Year-0)]*vpares$input$dat$waa[,as.character(w0_dat$Year-0)]))/1000
  
w0_dat = w0_dat %>% 
  mutate(Rec = rec, SSB = ssb,Biom=biom,Number=num,SSN = ssn)
# w0_dat$Year
# waa_dat$Year

plot(Weight~Number,data=w0_dat,log="x")
plot(Weight~SSB,data=w0_dat,log="x")
plot(Weight~Rec,data=w0_dat,log="x")
plot(Weight~Biom,data=w0_dat,log="x")
plot(Weight~SSN,data=w0_dat,log="x")


plot(Weight~Number,data=w0_dat,log="y")
plot(Weight~SSB,data=w0_dat,log="y")
plot(Weight~Rec,data=w0_dat,log="y")
plot(Weight~Biom,data=w0_dat,log="y")
plot(Weight~SSN,data=w0_dat,log="y")


plot(Rec~SSB,data=w0_dat,log="xy",col=Weight)

(g_tmp = ggplot(data=w0_dat,aes(x=SSB,y=Rec,colour=Weight))+
  geom_point(size=2)+
    scale_x_log10()+
    scale_y_log10()+
    scale_colour_gradient(high="darkorange",low="darkblue")
)

# cor(log(w0_dat$Rec),log(w0_dat$SSB_prev))
# lm(Weight~log(Number_prev),data=w0_dat) %>% summary
# lm(Weight~log(SSB_prev),data=w0_dat) %>% summary
# lm(Weight~log(Rec),data=w0_dat) %>% summary
# lm(Weight~log(Biom_prev),data=w0_dat) %>% summary


# full_w0 = glm(Weight~Number_prev+log(Number_prev),data=w0_dat,family=Gamma("identity"))
# full_w0 = glm(Weight~log(Number)+log(SSB)+log(Rec),data=w0_dat,family=Gamma("identity"))

# full_w0 = glm(Weight~Number+SSB+Rec,data=w0_dat,family=Gamma("identity"))
full_w0 = glm(Weight~Number+SSB*Rec+SSN*Rec,data=w0_dat,family=Gamma("log"))


# summary(full_w0)
# 
# dredge_w0 = dredge(full_w0, subset=!("log(Number)" &&"log(Rec)"))

dredge_w0 = dredge(full_w0, subset=!("Number" &&"Rec") & !("Number" &&"SSB") & !("Number" &&"SSN") & !("SSB" &&"SSN"))
head(dredge_w0,100) #470.9

# full_w0_2 = glm(Weight~Number+SSB+Rec,data=w0_dat,family=Gamma("identity"))

# dredge_w0_2 = dredge(full_w0_2, subset=!("Number" &&"Rec"))
# 
# 
# head(dredge_w0_2,100)

save(dredge_w0,file=savename("dredge_w0.rda"))

model_sel_w0 = model.sel(dredge_w0)
write.csv(as.data.frame(model_sel_w0),file=savename("AICc_table_w0.csv"))

mod_w0 = get.models(dredge_w0,subset=1)[[1]]
mod_w0$formula
# Weight ~ SSB + 1
summary(mod_w0)

# plot(predict(mod_w0)~w0_dat$SSB)
# 
# plot(predict(mod_w0),w0_dat$Weight)
# abline(0,1)

## figure weight at age 0 ----

# newdata_w0 = expand.grid(Number_prev=exp(seq(log(min(w0_dat$Number_prev)),log(max(w0_dat$Number_prev)),length=200))) %>%
#   as.data.frame()

newdata_w0 = expand.grid(
  # Rec=seq(min(w0_dat$Rec),max(w0_dat$Rec),length=200),
                         SSB=seq(0,max(w0_dat$SSB)*1.1,length=1000)) %>%
  as.data.frame()


tmp = predict(mod_w0,newdata=newdata_w0,type="response",se.fit=TRUE)
tmp$fit
newdata_w0 = newdata_w0 %>% mutate(Weight=tmp$fit,SE=tmp$se.fit) %>% 
  mutate(Upper=Weight+1.96*SE,Lower=Weight-1.96*SE)


# (g_w0 = ggplot(data=w0_dat,aes(x=SSB,y=Rec,colour=Weight))+
#     geom_contour_filled(data=newdata_w0,aes(z=Weight,fill = ..nlevel..),alpha=0.6)+
#     geom_point(size=2)+
#     scale_x_log10(expand = c(0.01, 0.01))+
#     scale_y_log10(expand = c(0.01, 0.01))+
#     scale_colour_gradient(high="darkorange",low="darkblue",name=bquote(italic(w[0])))+
#     scale_fill_gradient(high="darkorange",low="darkblue")+
#     guides(fill="none")+
#     ylab("Recruitment")+theme_bw(base_size=base_size)
# )


newdata_w0$Weight
(g_w0 = ggplot(data=w0_dat,aes(x=SSB,y=Weight))+
    geom_ribbon(data=newdata_w0,aes(ymax=Upper,ymin=Lower),alpha=0.5)+
    geom_path(data=newdata_w0,size=1)+
    geom_point(size=2)+
    # scale_x_log10(expand = c(0.01, 0.01))+
    # scale_y_log10(expand = c(0.01, 0.01))+
    # scale_colour_gradient(high="darkorange",low="darkblue",name=bquote(italic(w[0])))+
    # scale_fill_gradient(high="darkorange",low="darkblue")+
    # guides(fill="none")+
    ylab(bquote(italic(w[0])))+xlab("SSB (thousand ton)")+
    theme_bw(base_size=base_size)+
    scale_x_continuous(expand=c(0.02,0.02),limits=c(0,NA))+
    scale_y_continuous(expand=c(0.02,0.02),limits=c(0,NA))
)

# (g_w0 = ggplot(data=w0_dat,aes(x=SSB,y=Rec,z=Weight,colour=Weight))+
#     geom_contour(data=newdata_w0,aes(z=Weight),alpha=0.6)+
#     # geom_contour_filled(data=newdata_w0,aes(z=Weight,fill = ..level..),alpha=0.6)+
#     geom_point(size=2)+
#     scale_x_log10(expand = c(0.01, 0.01))+
#     scale_y_log10(expand = c(0.01, 0.01))+
#     scale_colour_gradient(high="darkorange",low="darkblue",name=bquote(italic(w[0])))+
#     # scale_fill_gradient(high="darkorange",low="darkblue")+
#     # guides(fill="none")+
#     ylab("Recruitment")+theme_bw(base_size=base_size)
# )

ggsave(g_w0,filename=savename("weight_age0.png"),dpi=600,height=100,width=150,unit="mm")


# (g_w0 = ggplot(data=NULL,aes(y=Weight,x=Number_prev))+
#     geom_ribbon(data=newdata_w0,aes(ymax=Upper,ymin=Lower),alpha=0.4)+
#     geom_point(data=w0_dat,size=point_size)+
#     # ylim(0,NA)+
#     geom_path(data=newdata_w0,size=path_size)+
#     # scale_colour_gradient(low="deepskyblue",high="sienna1",name="Abundance")+
#     theme_bw(base_size=base_size)+
#     ylab("Weight at age 0")+xlab("Abundance (billion)")+
#     scale_x_log10()+
#     expand_limits(y = 0)+
#     scale_y_continuous(expand = c(0.02, 0.02))
# )

save(mod_w_growth,file=savename("model_w_growth.rda"))
save(mod_w0,file=savename("model_w_age0.rda"))

save(g_w0,file=savename("weight_age0_graph.rda"))
ggsave(g_w0,filename=savename("weight_age0.png"),dpi=600,height=100,width=150,unit="mm")


## maturity modeling ----

# tranform [0,1] to (0,1)
# https://stats.stackexchange.com/questions/48028/beta-regression-of-proportion-data-including-1-and-0
vpares$input$dat$maa

delta_mat = function(Maturity,Maturity_prev) (Maturity-Maturity_prev)/(1-Maturity_prev)
y_trans = function(y0,N) (y0*(N-1)+0.5)/N
y0_trans = function(y,N) (y*N-0.5)/(N-1)
trans_mat = function(y0,Maturity_prev) Maturity_prev+y0*(1-Maturity_prev) 

# maa_dat = waa_dat2 %>% 
#   filter(Age>0 & Age<4) %>%
#   # filter(Year > min(Year)) %>%
#   mutate(y0 = delta_mat(Maturity,Maturity_prev)) %>%
#   mutate(N = n()) %>% 
#   mutate(y = y_trans(y0,N),Age=factor(Age)) 
# 
maa_dat = waa_dat2 %>% 
  filter(Age>0 & Age<4) %>%
  # filter(Year > min(Year)) %>%
  mutate(y0 = delta_mat(Maturity,Maturity_prev)) %>%
  mutate(N = n()) %>% 
  mutate(y = y_trans(y0,N),Age=factor(Age)) 

head(maa_dat)

waa_dat$Year %>% unique()
maa_dat$Year %>% unique()

plot(y~factor(Age),dat=maa_dat)

# vpares$naa
# maa_dat %>% filter(Age==3 & Maturity<1)
# 
# maa_dat3 = maa_dat %>% filter(Age==3) 


# full_mat = betareg(y~Age*Number_prev+Age*Cohort_prev+Age*Cohort_plus,
#                    data=maa_dat,link="logit",type="BC")
full_mat = betareg(y~Age*Number_prev+Age*Cohort_prev+Age*Cohort_plus,
                   data=maa_dat,link="cloglog",type="BC")

# glmmTMBに変更（0-1の応答変数を使えるため）

full_mat = glmmTMB(y0 ~ Age*Number_prev+Age*Cohort_prev+Age*Cohort_plus,
                   data=maa_dat,family=ordbeta)

?family_glmmTMB
?betareg
# full_mat = betareg(y~Age*DF_N+Age*DF_n+Age*DF_n_plus,
#                    data=maa_dat,link="logit",type="BC")

?betareg
# full_mat = betareg(y~Age*DF_N+Age*DF_n+Age*DF_n_plus,
#                    data=maa_dat,link="loglog",type="BC")
# loglog linkに変更(AICcが低くなる&logitだとfitが悪いように見えるため)

varying.link <- list(family = alist(logit = ordbeta("logit"),
                                    probit = ordbeta("probit"), cloglog = ordbeta("cloglog") ))

# dredge_mat = dredge(full_mat,
#                     subset=!("Number_prev" && "Cohort_prev") & 
#                       !("Cohort_plus" && "Cohort_prev") & !("Number_prev" && "Cohort_plus"),
#                     varying = varying.link)
# 
getAllTerms(full_mat)
dredge_mat = dredge(full_mat,
                    subset=!("cond(Number_prev)" && "cond(Cohort_prev)") & 
                      !("cond(Cohort_plus)" && "cond(Cohort_prev)") & !("cond(Number_prev)" && "cond(Cohort_plus)"),
                    varying = varying.link)

head(dredge_mat,20) 
# logit linkがよい
# nrow(maa_dat3)

mod_mat = get.models(dredge_mat,subset=1)[[1]]
mod_mat$call$formula
# y0 ~ Age + Cohort_plus + Age:Cohort_plus + 1
summary(mod_mat)
AICc(mod_mat)

# ?betareg
# 
# mod_mat_probit = update(mod_mat,link="probit")
# mod_mat_cloglog = update(mod_mat,link="cloglog")
# mod_mat_cauchit = update(mod_mat,link="cauchit")
# mod_mat_loglog = update(mod_mat,link="loglog")
# AICc(mod_mat,mod_mat_probit,mod_mat_cloglog,mod_mat_cauchit,mod_mat_loglog)
# 
# mod_mat <- mod_mat_loglog


model_sel_mat = model.sel(dredge_mat)
write.csv(as.data.frame(model_sel_mat),file=savename("AICc_table_maturity.csv"))

# mod_mat = get.models(dredge_mat,subset=1)[[1]]
# summary(mod_mat)
# AICc(mod_mat) #-608

# bind_cols(maa_dat,pred=predict(mod_mat)) %>% View

# predict(mod_mat)

### figure of maturity ----

newdata_mat = expand.grid(Age=unique(maa_dat$Age),N=unique(maa_dat$N),
                         Cohort_plus=seq(min(maa_dat$Cohort_plus),max(maa_dat$Cohort_plus),length=1000)) %>%
  as.data.frame()

# newdata_mat = expand.grid(Age=unique(maa_dat$Age),N=unique(maa_dat$N),
#                           DF_n_plus=seq(min(maa_dat$DF_n_plus),max(maa_dat$DF_n_plus),length=200)) %>%
#   as.data.frame()




?predict.betareg

# tmp = predict(mod_mat,newdata=newdata_mat)
# tmp
# 
# tmp2 = predict(mod_mat,newdata=newdata_mat,type="response")
# tmp2-tmp
# 

newdata_mat = newdata_mat %>% mutate(y=predict(mod_mat,newdata=newdata_mat,type="response")) %>%
  mutate(y0 = y0_trans(y,N)) %>%
  mutate(y0 = pmin(1,pmax(0,y0)))

# newdata_mat$y0%>% range
# newdata_mat$y0_l%>% range

# newdata_mat$y %>% range
# newdata_mat$y0 %>% range

# temp2 = newdata_mat %>% filter(Age=="1")
# plot(temp2$Number_prev,temp2$Maturity)

# tmp = predict(mod_mat)
# tmp %>% summary
# temp = bind_cols(maa_dat,data.frame(predd=tmp)) %>% filter(Age=="1")

# plot(temp$Number_prev,temp$predd)




(g_mat = ggplot(data=NULL,aes(x=Cohort_plus))+
    geom_point(data=maa_dat,aes(y=y0,colour=Age),size=point_size)+
    # xlim(0,NA)+ylim(0,NA)+
    geom_path(data=newdata_mat,aes(y=y0,colour=Age,group=Age),size=path_size)+
    scale_colour_brewer(palette="Dark2")+
    theme_bw(base_size=base_size)+
    # scale_x_log10(expand = c(0.02, 0.02))+
    xlab(bquote(italic(n)*'+'))+ylab(bquote(Delta*italic(g)))+
    expand_limits(y = 0)+
    # scale_x_log10(expand = c(0.02, 0.02)) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0.02, 0.02))
)

ggsave(g_mat,filename=savename("maturity_growth.png"),dpi=600,height=100,width=150,unit="mm")
save(g_mat,file=savename("maturation_graph.rda"))


## generate combined figure ----

g1 = gridExtra::grid.arrange(g_w0 + ggtitle("(a)"),
                             g_wg + ggtitle("(b)"),
                             g_mat + ggtitle("(c)"),nrow=3)

ggsave(g1,filename=savename("individual_dens_effect.png"),dpi=600,unit="mm",height=240,width=150)
save(g1,file=savename("individual_dens_effect.rda"))


g1_wide = gridExtra::grid.arrange(g_w0 + ggtitle("(a)"),
                             g_wg + ggtitle("(b)"),
                             g_mat + ggtitle("(c)"),nrow=1)

ggsave(g1_wide,filename=savename("individual_dens_effect_wide.png"),dpi=600,unit="mm",height=80,width=320)
save(g1_wide,file=savename("individual_dens_effect_wide.rda"))

## analyze density-independent model ----

summary(mod_w0)  #SSBの方が影響大きい
mod_w0_di = update(mod_w0, formula=~1)
summary(mod_w0_di)

summary(mod_w_growth)
mod_w_growth_di = update(mod_w_growth, formula=~Weight_prev)
summary(mod_w_growth_di)

summary(mod_mat)
mod_mat_di = update(mod_mat, formula=~Age)
summary(mod_mat_di)

AICc(mod_mat,mod_mat_di)
AICc(mod_w0,mod_w0_di)
AICc(mod_w_growth,mod_w_growth_di)


###  SR relationship ----

SRdata=get.SRdata(vpares,years=unique(waa_dat$Year))
SRdata$SSB <- SRdata$SSB/1000
SRdata$R <- SRdata$R/1000

# resHS = fit.SR(SRdata,SR="HS",AR=0,out.AR=FALSE)
# resHS$pars

resBH = fit.SR(SRdata,SR="BH",AR=0,out.AR=FALSE)
resBH$pars

# resBH0 = fit.SR(SRdata,SR="BH",AR=0,out.AR=FALSE)
# resBH0$pars
# 
# resRI = fit.SR(SRdata,SR="RI",AR=1,out.AR=FALSE)
# resRI$pars
# resRI$pred


(g_BH = SRplot_gg(resBH))
# (g_HS = SRplot_gg(resHS))
# (g_RI = SRplot_gg(resRI))

# 1/resHS$pars$a

# predict(resBH$opt)

# c(resHS$AICc,resBH$AICc)

# using nls() CIを求めるため


# bh_log = function(loga,logb,x) log(exp(loga)*x)-log(1+exp(logb)*x)
# 
# hs = function(loga,transb,x) {
#   a = exp(loga)
#   b = min(SRdata$SSB) + (max(SRdata$SSB)-min(SRdata$SSB))/(1+exp(-transb))
#   return( a*min(b,x) )
# }

# init = list(loga=resBH$opt$par[1],logb=resBH$opt$par[2])
# bh_res = nls(log(R/SSB)~log(exp(loga))-log(1+exp(logb)*SSB),start=init,data=SRdata)

# init = list(loga=resHS$opt$par[1],transb=resHS$opt$par[2])
# hs_res = nls(log(R)~log(exp(loga))-log(1+exp(logb)*SSB),start=init,data=SRdata)
# 
# ?nls

# pred_HS = resHS$pred %>% as.data.frame() 
pred_BH = resBH$pred %>% as.data.frame() 
# tmp = predFit(bh_res,newdata=pred_BH,se.fit=TRUE, interval = "confidence", level= 0.95, adjust="none",k=100)
# tmp$se.fit
# tmp$fit
# plot(tmp$fit[,3])
# 
# pred_BH = pred_BH %>% bind_cols(as.data.frame(tmp$fit)) %>% 
#   mutate(pred = exp(fit)*SSB,upper=exp(upr)*SSB,lower=exp(lwr)*SSB)

# (g_HS = ggplot(data=NULL,aes(x=SSB))+
#     # geom_ribbon(data=pred_BH,aes(ymax=upper,ymin=lower),alpha=0.4)+
#     geom_path(data=pred_HS,aes(y=R),size=1)+
#     geom_point(data=SRdata,aes(y=R),size=2)+
#     theme_bw(base_size=12)+
#     xlab("SSB (thousand ton)")+ylab("Recruits (billion)") +
#     ggtitle("(a)")
# )
# 
# ggsave(g_HS,filename=savename("HS.png"),dpi=600,unit="mm",height=100,width=150)

(g_BH = ggplot(data=NULL,aes(x=SSB))+
    # geom_ribbon(data=pred_BH,aes(ymax=upper,ymin=lower),alpha=0.4)+
    geom_path(data=pred_BH,aes(y=R),size=1)+
    geom_point(data=SRdata,aes(y=R),size=2)+
    theme_bw(base_size=12)+
    xlab("SSB (thousand ton)")+ylab("Recruits (billion)") +
    ggtitle("(a)")
)

ggsave(g_BH,filename=savename("BH.png"),dpi=600,unit="mm",height=100,width=150)


# Scenario A: use five years average ----

# HS = function(x,a=as.numeric(resHS$pars["a"]),b=as.numeric(resHS$pars["b"])) SRF_HS(x,a,b)
# HS(100)
# HS2(100)
# SRF_HS(100,)

BH = function(x,a=as.numeric(resBH$pars["a"]),b=as.numeric(resBH$pars["b"])) SRF_BH(x,a,b)
BH(100)

vpares$input$dat$waa

fc_year = as.character(2017:2021)

Fcurrent = rowMeans(vpares$faa[,fc_year])
maa_A = rowMeans(vpares$input$dat$maa[,fc_year])
M_A = rowMeans(vpares$input$dat$M[,fc_year])
waa_A = rowMeans(vpares$input$dat$waa[,fc_year])

# resB0_A = frasyr::calc_steepness(SR="HS",rec_pars=resHS$pars,
#                                  M=M_A,maa=maa_A,waa=waa_A,faa=Fcurrent)

resB0_A = frasyr::calc_steepness(SR="BH",rec_pars=resBH$pars,
                                 M=M_A,maa=maa_A,waa=waa_A,faa=Fcurrent)

resB0_A


# Scenario B: use all years average ----

maa_B = rowMeans(vpares$input$dat$maa)
M_B = rowMeans(vpares$input$dat$M)
waa_B = rowMeans(vpares$input$dat$waa)

# resB0_B = frasyr::calc_steepness(SR="HS",rec_pars=resHS$pars,
#                                  M=M_B,maa=maa_B,waa=waa_B,faa=Fcurrent)

resB0_B = frasyr::calc_steepness(SR="BH",rec_pars=resBH$pars,
                                 M=M_B,maa=maa_B,waa=waa_B,faa=Fcurrent)

res_B0 = full_join(resB0_A %>% mutate(Scenario="A"),
                   resB0_B %>% mutate(Scenario="B"))
res_B0



# Scenario C: density-independent model ----

summary(mod_w0_di)
waa_C = mean(predict(mod_w0_di,type="response")) # age 0 weight

summary(mod_w_growth_di)

for (i in 1:A) {
  tmp = predict(mod_w_growth_di,newdata= data.frame(Weight_prev=waa_C[i]))
  waa_C = c(waa_C,tmp)
}

cbind(waa_A,waa_B,waa_C)

tmp = maa_dat %>% mutate(pred=predict(mod_mat_di)) %>%
  mutate(pred=y0_trans(pred,N=nrow(maa_dat))) %>% 
  group_by(Age) %>%
  summarise(value = mean(pred)) %>%
  mutate(age = as.numeric(Age))
pred_mat <- 0
for (i in 1:nrow(tmp)) {
  pred_mat <- c(pred_mat,trans_mat(tmp$value[i], rev(pred_mat)[1]))
}
# pred_mat <- sapply(1:nrow(tmp), function(i) trans_mat(tmp$value[i], c(0,tmp$value[-3])[i]))
# 
# trans_mat(tmp$value[3],pred_mat[2])

maa_C <- c(pred_mat,rep(1,3))

# maa_C <- maa_B
# maa_C[names(maa_C) %in% tmp$Age] <- pred_mat
# 
cbind(maa_A,maa_B,maa_C)
M_C <- M_B

# resB0_C = frasyr::calc_steepness(SR="HS",rec_pars=resHS$pars,
#                                  M=M_C,maa=maa_C,waa=waa_C,faa=Fcurrent)

resB0_C = frasyr::calc_steepness(SR="BH",rec_pars=resBH$pars,
                                 M=M_C,maa=maa_C,waa=waa_C,faa=Fcurrent)

res_B0 = res_B0 %>% full_join(.,
                   resB0_C %>% mutate(Scenario="C"))
res_B0


# Scenario D: density-dependent model ----

source("~/git/masaba_dd/source.R")

## Step1: Calculating equilibrium ----

# ここから先はやっていない 2023/09/21)

# 
# rc = seq(0.001,max(res_B0$R0),length=200)
# M = 0.4
# x = 1
# x = resB0_C$Fmsy2F
# 
# faa = x*Fcurrent
# resSR = resBH
# 
# model_mat = mod_mat
# model_w0 = mod_w0
# model_wg = mod_w_growth

# model_mat = mod_mat_di
# model_w0 = mod_w0_di
# model_wg = mod_w_growth_di

x = seq(0.00,10.00,by=0.1)
EQdata <- AGEdata <- data.frame()
# eqres$EQdata
# eqres$AGEdata

# undebug(calc_eq)
# undebug(calc_SPR)
# undebug(opt_w0)
# undebug(pred_mat)
# 
i=1
for (i in 1:length(x)) {
  faa = Fcurrent*x[i]
  eqres = calc_eq(faa=faa,model_mat = mod_mat,model_w0 = mod_w0,model_wg = mod_w_growth,resSR=resBH)
  EQdata = EQdata %>% bind_rows(eqres$EQdata %>% mutate(Fmulti = x[i]))
  AGEdata = AGEdata %>% bind_rows(eqres$AGEdata %>% mutate(Fmulti = x[i]))
  print(paste0(i, ", Fmulti=",round(x[i],2),", SY=",round(eqres$EQdata$SY), ", SSB=",round(eqres$EQdata$SSB)))
}

calc_eq

plot(EQdata$SSB/EQdata$SPR0~EQdata$SSB)

write.csv(EQdata,file=savename("EQdata.csv"))
write.csv(AGEdata,file=savename("AGEdata.csv"))

# AGEdata = read.csv(file=savename("EQdata.csv"))
# EQdata = AGEdata %>% group_by(Fmulti) %>% 
#   summarise(B=sum(baa),SSB=sum(ssb),SY=sum(catch),
#             SPR0 = ref.F(Fcurrent=faa,waa=waa,maa=maa,M=M,waa.catch=waa,Pope=TRUE,F.range=1)$spr0
#               )
#   
# EQdata = EQdata %>% full_join(
#   AGEdata %>% filter(age==0) %>% 
#   rename(R = naa) %>% 
#   dplyr::select(Fmulti,R)
# )
# 
# EQdata = EQdata %>% mutate(SPR=SSB/R) %>%
#   mutate(pSPR = SPR/SPR0)

write.csv(EQdata,file=savename("EQdata.csv"))


## Step2: Estimating MSY ----

xl = x[which.max(EQdata$SY)-1]
xu = x[which.max(EQdata$SY)+1]

obj_msy = function(x,out=FALSE) {
  faa = Fcurrent*x
  eqres = calc_eq(faa=faa,model_mat = mod_mat,model_w0 = mod_w0,model_wg = mod_w_growth,resSR=resBH)
  if (out==FALSE) return( -eqres$EQdata$SY )
  if (out==TRUE) return (eqres)
}

opt_msy = optimize(obj_msy,c(xl,xu))
temp = obj_msy(opt_msy$minimum,out=TRUE)
resB0_D = temp$EQdata %>% 
  dplyr::select(-YPR) %>%
  rename(Bmsy=B,SBmsy=SSB,Rmsy=R,MSY=SY,SPRmsy=SPR) %>%
  mutate(Scenario="D",Fmsy2F=opt_msy$minimum)

resB0_D$SBmsy/resB0_D$Rmsy

resB0_D = resB0_D %>% bind_cols(
  filter(EQdata,Fmulti==0) %>% 
  # dplyr::select(-YPR) %>%
  rename(B0=B,SB0=SSB,R0=R) %>%
    dplyr::select(B0,SB0,R0)
) %>% mutate(h = 1-SBmsy/SB0)

resB0_msy = full_join(res_B0,resB0_D)
resB0_msy = resB0_msy %>% mutate(pSPR = SPRmsy/SPR0) %>%
  dplyr::select(Scenario,everything())

write.csv(resB0_msy,file=savename("resB0_msy.csv"))


### figure maturity at age and weight at age ----

out_D = temp$AGEdata
out_D$maa
mgdat = data.frame(age=0:A,value=maa_A,type="Maturity",scenario="A") %>%
  full_join(data.frame(age=0:A,value=maa_B,type="Maturity",scenario="B")) %>%
  full_join(data.frame(age=0:A,value=maa_C,type="Maturity",scenario="C")) %>%
  full_join(data.frame(age=0:A,value=out_D$maa,type="Maturity",scenario="D")) %>% 
  full_join(data.frame(age=0:A,value=waa_A,type="Weight",scenario="A")) %>%
  full_join(data.frame(age=0:A,value=waa_B,type="Weight",scenario="B")) %>%
  full_join(data.frame(age=0:A,value=waa_C,type="Weight",scenario="C")) %>%
  full_join(data.frame(age=0:A,value=out_D$waa,type="Weight",scenario="D"))

(g_mg = ggplot(data=mgdat,aes(x=age,y=value,colour=scenario,linetype=scenario))+
  geom_path(size=1)+
  facet_wrap(vars(type),scales="free_y")+
    theme_bw(base_size=base_size)+
    theme(legend.position="top")+
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),name="Scenario")+
    scale_colour_brewer(palette="Dark2",name="Scenario")+ylim(0,NA)+
    xlab("Age")+ylab("")
)

ggsave(g_mg,filename=savename("maturity-weight.png"),dpi=600,unit="mm",height=90,width=180)


### Stock-recruitment relationship and replacement line ----  

replace_dat = EQdata %>% 
  mutate(R2 = SSB/SPR0) %>%
  filter(SSB < max(pred_BH$SSB))
  

sprres = get.SPR(vpares)
ysdata = sprres$ysdata %>% mutate(Year = as.numeric(rownames(sprres$ysdata))) %>%
  mutate(SSB = as.numeric(colSums(vpares$ssb))/1000) %>% 
  filter(Year < max(Year)) %>% 
  mutate(R2 = SSB/SPR0)


# EQdata$SSB %>% range

replace_dat_all2 = EQdata %>%
  dplyr::select(SSB,SPR0,SY) %>% 
  mutate(Scenario="D") %>%
    mutate(R2 = SSB/SPR0)


replace_dat_all = EQdata %>% 
  filter(Fmulti > 0) %>% 
  dplyr::select(SSB,SPR0,SY) %>% 
  mutate(Scenario="D")

ssbc = seq(0,max(replace_dat_all$SSB),length=1000)

replace_dat_all =replace_dat_all %>%
  full_join(data.frame(SSB = ssbc,SPR0 = resB0_A$SPR0,Scenario ="A") ) %>%
  full_join(data.frame(SSB = ssbc,SPR0 = resB0_B$SPR0,Scenario ="B") ) %>%
  full_join(data.frame(SSB = ssbc,SPR0 = resB0_C$SPR0,Scenario ="C") ) %>%
  arrange(Scenario,SSB)


replace_dat_all3 = replace_dat_all

replace_dat_all = replace_dat_all %>%
  mutate(R2 = SSB/SPR0) %>% 
  filter(SSB < max(pred_BH$SSB))



(g_replace = ggplot(data=NULL,aes(x=SSB))+
    geom_path(data=replace_dat_all,aes(y=R2,colour=Scenario,linetype=Scenario),size=path_size)+
    geom_point(data=ysdata,aes(y=R2),size=point_size)+
    scale_colour_brewer(palette="Dark2")+theme_bw(base_size=base_size)+
    theme(legend.position="top")+
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("SSB (thousand ton)")+ylab("Recruits (billion)")
)

# ggsave(g_replace,filename=savename("replacement.png"),dpi=600,unit="mm",height=100,width=150)

library(RColorBrewer)

(g_BH2 = ggplot(data=NULL,aes(x=SSB))+
    # geom_ribbon(data=pred_BH,aes(ymax=upper,ymin=lower),alpha=0.4)+
    geom_path(data=pred_BH,aes(y=R),size=path_size)+
    geom_point(data=SRdata,aes(y=R),size=point_size)+
    geom_path(data=replace_dat,aes(y=R2),linetype="dashed",colour=brewer.pal(4, "Dark2")[4],size=1)+
    geom_vline(aes(xintercept=resB0_D$SBmsy),linetype="dotted",colour=brewer.pal(4, "Dark2")[4],size=1)+
    theme_bw(base_size=12)+
    xlab("SSB (thousand ton)")+ylab("Recruits (billion)")
)

ggsave(g_BH2,filename=savename("BH_replace.png"),dpi=600,unit="mm",height=100,width=150)

(g_replace2 = ggplot(data=NULL,aes(x=SSB))+
    geom_path(data=replace_dat_all,aes(y=R2,colour=Scenario,linetype=Scenario),size=path_size)+
    geom_point(data=ysdata,aes(y=R2),size=point_size)+
    scale_colour_brewer(palette="Dark2")+theme_bw(base_size=base_size)+
    # theme(legend.position="top")+
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("SSB (thousand ton)")+ylab("Recruits (billion)")+
    theme(legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.key.size=unit(1, 'cm'),
          legend.position=c(0.02,0.98),legend.justification = c(0,1))
)

ggsave(g_replace2,filename=savename("replacement.png"),dpi=600,unit="mm",height=100,width=150)

g_BH_replace = gridExtra::grid.arrange(g_BH2+ggtitle("(a)"),g_replace2+ggtitle("(b)"),nrow=1)

ggsave(g_BH_replace,filename=savename("BH&replace.png"),dpi=600,unit="mm",height=90,width=240)


## spawners per spawner ---

replace_dat_all$SSB %>% range

sydata = replace_dat_all3 %>% 
  mutate(R = BH(SSB)) %>%
  mutate(RS = R*SPR0) %>%
  mutate(SPS = RS/SSB) %>%
  filter(SSB>0)


# sydata %>% filter(Scenario=="D") %>% View

(g_sps = ggplot(sydata,aes(x=SSB,y=SPS,colour=Scenario,linetype=Scenario))+
  geom_path(size=path_size)+    
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("SSB (thousand ton)")+ylab("Spawners per spawner")+
    theme_bw(base_size=base_size)+
  theme(legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=9),
          legend.key.size=unit(1, 'cm'),
          legend.position=c(0.98,0.98),legend.justification = c(1,1))
)

ggsave(g_sps,filename=savename("spawners_per_spawner.png"),dpi=600,unit="mm",height=100,width=150)


(g_sps2 = ggplot(sydata,aes(x=SSB,y=SPS,colour=Scenario,linetype=Scenario))+
    geom_path(size=path_size)+    
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(nrow=1))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(nrow=1))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("SSB (thousand ton)")+ylab("Spawners per spawner")+
    theme_bw(base_size=base_size)+
    theme(legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=9),
          legend.key.size=unit(1, 'cm'),
          legend.position=c(0.98,0.98),legend.justification = c(1,1))
)

ggsave(g_sps2,filename=savename("spawners_per_spawner2.png"),dpi=600,unit="mm",height=100,width=150)


## Sustainable yield  ----

### Scenario C ----

# EQdata_all = EQdata %>% mutate(Scenario="D")
# AGEdata_all = AGEdata %>% mutate(Scenario="D")

EQdata_C <- AGEdata_C <- data.frame()

for (i in 1:length(x)) {
  faa = Fcurrent*x[i]
  eqres = calc_eq_fix(faa=faa,waa=waa_C,maa=maa_C,resSR=resBH)
  EQdata_C = EQdata_C %>% bind_rows(eqres$EQdata %>% mutate(Fmulti = x[i]))
  AGEdata_C = AGEdata_C %>% bind_rows(eqres$AGEdata %>% mutate(Fmulti = x[i]))
  print(paste0(i, ", Fmulti=",round(x[i],2),", SY=",round(eqres$EQdata$SY), ", SSB=",round(eqres$EQdata$SSB)))
}

### Scenario B ----

# EQdata_all = EQdata %>% mutate(Scenario="D")
# AGEdata_all = AGEdata %>% mutate(Scenario="D")

EQdata_B <- AGEdata_B <- data.frame()

for (i in 1:length(x)) {
  faa = Fcurrent*x[i]
  eqres = calc_eq_fix(faa=faa,waa=waa_B,maa=maa_B,resSR=resBH)
  EQdata_B = EQdata_B %>% bind_rows(eqres$EQdata %>% mutate(Fmulti = x[i]))
  AGEdata_B = AGEdata_B %>% bind_rows(eqres$AGEdata %>% mutate(Fmulti = x[i]))
  print(paste0(i, ", Fmulti=",round(x[i],2),", SY=",round(eqres$EQdata$SY), ", SSB=",round(eqres$EQdata$SSB)))
}

### Scenario A ----

EQdata_A <- AGEdata_A <- data.frame()

x2 = seq(0,2,length=101)

for (i in 1:length(x2)) {
  faa = Fcurrent*x2[i]
  eqres = calc_eq_fix(faa=faa,waa=waa_A,maa=maa_A,resSR=resBH)
  EQdata_A = EQdata_A %>% bind_rows(eqres$EQdata %>% mutate(Fmulti = x2[i]))
  AGEdata_A = AGEdata_A %>% bind_rows(eqres$AGEdata %>% mutate(Fmulti = x2[i]))
  print(paste0(i, ", Fmulti=",round(x2[i],2),", SY=",round(eqres$EQdata$SY), ", SSB=",round(eqres$EQdata$SSB)))
}


### combined ----

EQdata_all = EQdata_A %>% mutate(Scenario="A") %>% 
  full_join(EQdata_B %>% mutate(Scenario="B")) %>% 
  full_join(EQdata_C %>% mutate(Scenario="C")) %>% 
  full_join(EQdata %>% mutate(Scenario="D"))

write.csv(EQdata_all,file=savename("EQdata_all.csv"))
  
AGEdata_all = AGEdata_A %>% mutate(Scenario="A") %>% 
  full_join(AGEdata_B %>% mutate(Scenario="B")) %>% 
  full_join(AGEdata_C %>% mutate(Scenario="C")) %>% 
  full_join(AGEdata %>% mutate(Scenario="D"))

write.csv(AGEdata_all,file=savename("AGEdata_all.csv"))

obs_sy = data.frame(Year=as.numeric(names(colSums(vpares$ssb))),
                    SSB=as.numeric(colSums(vpares$ssb))/1000,
                    SY = as.numeric(colSums(vpares$input$dat$caa*vpares$input$dat$waa))/1000) %>%
  filter(Year<max(Year))

obs_range = expand.grid(SSB=range(SRdata$SSB),SY=c(0,max(EQdata_all$SY)))[c(1,2,4,3),]


(g_sy = ggplot(EQdata_all,aes(x=SSB,y=SY))+
    geom_polygon(data=obs_range,alpha=0.3)+
    # geom_point(data=obs_sy,size=point_size)+
    geom_path(aes(colour=Scenario,linetype=Scenario),size=path_size)+
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("SSB (thousand ton)")+ylab("Sustainable yield (thousand ton)")+
    theme_bw(base_size=base_size)+
    theme(legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=9),
          legend.key.size=unit(1, 'cm'),
          legend.position=c(0.98,0.98),legend.justification = c(1,1))
)

ggsave(g_sy,filename=savename("sustainable_yield.png"),dpi=600,unit="mm",height=100,width=150)


(g_sy2 = ggplot(EQdata_all,aes(x=SSB,y=SY))+
    geom_polygon(data=obs_range,alpha=0.3)+
    # geom_point(data=obs_sy,size=point_size)+
    geom_path(aes(colour=Scenario,linetype=Scenario),size=path_size)+
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("SSB (thousand ton)")+ylab("Sustainable yield (thousand ton)")+
    theme_bw(base_size=base_size)+
    theme(legend.position="none")
)

ggsave(g_sy2,filename=savename("sustainable_yield2.png"),dpi=600,unit="mm",height=100,width=150)


EQdata_all2 = EQdata_all %>% left_join(resB0_msy,by="Scenario") %>%
  mutate(relSY = SY/MSY,relSSB=SSB/SB0)

(g_fsy = ggplot(EQdata_all2,aes(x=Fmulti,y=SY))+
    geom_path(aes(colour=Scenario,linetype=Scenario),size=path_size)+    
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("F relative to Fcurrent")+ylab("Sustainable yield (thousand ton)")+
    theme_bw(base_size=base_size)+
    theme(legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=9),
          legend.key.size=unit(1, 'cm'),
          legend.position=c(0.98,0.98),legend.justification = c(1,1))
)

ggsave(g_fsy,filename=savename("sustainable_yield2F.png"),dpi=600,unit="mm",height=100,width=150)


(g_fsy2 = ggplot(EQdata_all2,aes(x=Fmulti,y=SY))+
    geom_path(aes(colour=Scenario,linetype=Scenario),size=path_size)+    
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("F relative to Fcurrent")+ylab("Sustainable yield (thousand ton)")+
    theme_bw(base_size=base_size)+
    theme(legend.position="none")
)

ggsave(g_fsy2,filename=savename("sustainable_yield2F2.png"),dpi=600,unit="mm",height=100,width=150)


(g_relSSB = ggplot(EQdata_all2,aes(x=Fmulti,y=relSSB))+
    geom_path(aes(colour=Scenario,linetype=Scenario),size=path_size)+    
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("F relative to Fcurrent")+ylab("SSB relative to SSB0")+
    theme_bw(base_size=base_size)+
    theme(legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=9),
          legend.key.size=unit(1, 'cm'),
          legend.position=c(0.98,0.98),legend.justification = c(1,1))
)

ggsave(g_relSSB,filename=savename("relative_SSB.png"),dpi=600,unit="mm",height=100,width=150)

(g_relSSB2 = ggplot(EQdata_all2,aes(x=Fmulti,y=relSSB))+
    geom_path(aes(colour=Scenario,linetype=Scenario),size=path_size)+    
    scale_linetype_manual(values=c("dashed","dotted","dotdash","solid"),guide=guide_legend(ncol=2))+
    scale_colour_brewer(palette="Dark2",guide=guide_legend(ncol=2))+
    scale_x_continuous(expand=c(0.02,0.02))+
    scale_y_continuous(expand=c(0.02,0.02))+
    xlab("F relative to Fcurrent")+ylab("SSB relative to SSB0")+
    theme_bw(base_size=base_size)+
    theme(legend.position="none")
)

ggsave(g_relSSB2,filename=savename("relative_SSB2.png"),dpi=600,unit="mm",height=100,width=150)


g_eq = gridExtra::grid.arrange(g_sps2+ggtitle("(a)"),
                               g_sy2+ggtitle("(b)"),
                               g_fsy2+ggtitle("(c)"),
                               g_relSSB2+ggtitle("(d)"),nrow=2)


ggsave(g_eq,filename=savename("SPS&SustainableYield.png"),dpi=600,unit="mm",height=180,width=240)


head(replace_dat_all)
# sort(colSums(vpares$ssb))

# waa = waa_A
# maa = maa_A

obj_fun_di = function(x,s,waa,maa,SR="BH",out=FALSE) {
  
  if (SR=="BH") r = BH(s) else  r = HS(s)

  tmp = calc_rel_naa(faa=Fcurrent*x)
  naa = r*tmp
  # N = sum(naa)

  ssb_x = sum(naa*waa*maa)
  diff = ((s-ssb_x)/s)^2
  
  if (out==FALSE) return(diff) 
  if (out==TRUE) return(list(naa=naa,waa=waa,maa=maa,faa=Fcurrent*x)) 
}


ssb = seq(0.1, max(res_B0$SB0),length=200)

eqdat <- data.frame()

# j <- 1
# i = 1
# 
# obj_fun_di(0.1,0.1,waa=WAA,maa=MAA)
# # 
# sapply(seq(0,2,by=0.1), function(x) obj_fun_di(x,s=0.1,waa=WAA,maa=MAA))

j = 4

for(j in 1:4) {
  for (i in 1:length(ssb)) {
    print(i)
    
    if (j < 4) {
      WAA = case_when(j==1 ~ waa_A,
                      j==2 ~ waa_B,
                      TRUE ~ waa_C)
      MAA = case_when(j==1 ~ maa_A,
                      j==2 ~ maa_B,
                      TRUE ~ maa_C)
      opt = optimize(obj_fun_di,c(0,100),s=ssb[i],waa=WAA,maa=MAA)
      # as.data.frame(opt)
    } else {
      opt = optimize(obj_fun,c(0,100),s=ssb[i])
    }
    Fmulti = ifelse(opt$objective < 0.001 & opt$minimum>0.001,opt$minimum,0)
    
    if (j < 4) {
      out = obj_fun_di(x=Fmulti,s=ssb[i],waa=WAA,maa=MAA,out=TRUE)
    } else {
      out = obj_fun(x=Fmulti,s=ssb[i],out=TRUE)
    }
    out = as.data.frame(out)
    
    res_refF = frasyr::ref.F(Fcurrent=out$faa,Pope=TRUE,max.age=100,
                             maa=out$maa,waa=out$waa,M=rep(0.4,A+1),
                             waa.catch=out$waa,min.age=0,plot=FALSE,pSPR=0,F.range=1)
    
    yprspr = res_refF$ypr.spr[1,c("ypr","pspr")]
    SY = as.numeric(yprspr[1]*out$naa[1])
    SPR0 =as.numeric(res_refF$spr0)
    
    Biomass = sum(out$naa*out$waa)
    SSB = sum(out$naa*out$maa*out$waa)
    Rec = out$naa[1]
    
    naa0 = calc_rel_naa()
    SPS = sum(Rec*naa0*out$waa*out$maa)/SSB
    
    eqdat = eqdat %>% bind_rows(., 
                                data.frame(Fmulti=Fmulti,Biomass=Biomass,SSB=SSB,Rec=Rec,SPR0=SPR0,SPS=SPS,Yield=SY) %>%
                                  mutate(Scenario=c("A","B","C","D")[j])
    )
    # if (Fmulti==0) break
  }
}

eqdat$SSB

write.csv(eqdat,file=savename("eqdat.csv"))

res_B0msy

(g_eq1 = ggplot(data=eqdat,aes(x=Fmulti,y=Yield,colour=Scenario,linetype=Scenario))+
  geom_path()
)

eqdat %>% filter(Scenario=="A")
eqdat %>% filter(Scenario=="B")
eqdat %>% filter(Scenario=="D")

(g_eq2 = ggplot(data=eqdat,aes(x=Fmulti,y=SSB,colour=Scenario,linetype=Scenario))+
    geom_path()
)

(g_eq3 = ggplot(data=eqdat,aes(x=SSB,y=SPS,colour=Scenario,linetype=Scenario))+
    geom_path()+
    scale_x_continuous(limits=c(0,2500))
)

(g_eq4 = ggplot(data=eqdat,aes(x=SSB,y=Yield,colour=Scenario,linetype=Scenario))+
    geom_path()
)

View(eqdat)

i = 10

ssb[i]


plot(dd_eq$SSB,dd_eq$Rec)

res_B0

dd_eq %>% filter(Yield==max(Yield))

res_B0

summary(mod_mat)

### draw replacement line ----

# range(waa_dat$Number_prev,na.rm = TRUE)
# 
# model_w0 = mod_w0_di
# model_wg = mod_w_growth_di
# model_mat = mod_mat_di

# model_w0 = mod_w0
# model_wg = mod_w_growth
# model_mat = mod_mat

# summary(model_mat)

# dd = TRUE



calc_replace = function(dd) {
  
  if (isTRUE(dd)) {
    model_w0 = mod_w0
    model_wg = mod_w_growth
    model_mat = mod_mat
  } else {
    model_w0 = mod_w0_di
    model_wg = mod_w_growth_di
    model_mat = mod_mat_di
  }
  
  Rc = seq(min(vpares$naa[1,]/1000)*0.001,0.5*max(as.numeric(vpares$naa[1,])/1000),length=200)
  Fcurrent = rowMeans(vpares$faa[,as.character(2016:2020)])
 
  AGEdata = data.frame()
  SPRdata = data.frame()
  # l=50
  for (l in 1:length(Rc)) {
    r = Rc[l] # recruitment
    res = calc_SPR(rec=r,model_w0,model_wg,model_mat,faa=Fcurrent*0,M=0.4)
    AGEdata = bind_rows(AGEdata,res$AGEdata %>% mutate(ID=l))
    SPRdata = bind_rows(SPRdata,res$SPRdata %>% mutate(ID=l))
  }
  return(list(AGEdata=AGEdata,SPRdata=SPRdata))
}


# model=mod_mat_di
# naa=res$AGEdata$naa
# waa=res$AGEdata$waa
# N=nrow(model$model)

replace_di = calc_replace(dd=FALSE)
replace_dd = calc_replace(dd=TRUE)

replace_di$SPRdata$SPR %>% unique()
replace_dd$SPRdata$SPR 

replace_dat = replace_di$SPRdata %>% mutate(type="Independent") %>%
  full_join(replace_dd$SPRdata %>% mutate(type="Dependent"))

(g_replacement = ggplot(filter(replace_dat,SSB<2000),aes(x=SSB,y=R,colour=type,group=type))+
  geom_path()
)


obj_msy = function(x,Fcurrent,dd=FALSE,resSR,out=FALSE) {
  if (isTRUE(dd)) {
    model_w0 = mod_w0
    model_wg = mod_w_growth
    model_mat = mod_mat
  } else {
    model_w0 = mod_w0_di
    model_wg = mod_w_growth_di
    model_mat = mod_mat_di
  }
  tmp = calc_eq(faa=Fcurrent*x,model_w0,model_wg,model_mat,resSR=resSR)
  if (out==FALSE) {
    return( -tmp$EQdata$SY[1] )
  } else {
    return( tmp )
  }
}

x = 2
faa = Fcurrent*x
resSR = resBH
dd = TRUE
M = 0.4

res_B0msy

obj_msy(1, Fcurrent=Fcurrent,dd=TRUE,resSR=resBH)

det_init = sapply(seq(0.1,10,length=100), function(x) obj_msy(x, Fcurrent=Fcurrent,dd=FALSE,resSR=resBH))
det_init



msy_dd = optimize(obj_msy,interval=c(0.01,10),Fcurrent=Fcurrent,dd=TRUE,resSR=resBH)
msy_dd
temp = obj_msy(msy_dd$minimum,Fcurrent=Fcurrent,dd=TRUE,resSR=resBH,out=TRUE)
temp$EQdata
temp$EQdata
shape = MASS::gamma.shape(mod_w0_di)
shape$alpha
mean=as.numeric(predict(mod_w0_di)[1])
tmp = rgamma(100000,shape=shape$alpha,scale=mean/shape$alpha)
hist(tmp)
summary(tmp)

hist(1/tmp)
summary(1/tmp)

mean(1/tmp)/(1/mean)
summary(mod_w0_di)

resBH$pred
resHS$pred
g_BH = SRplot_gg(resBH)

replace_di$SSB/replace_di$R
cbind(replace_dd$SSB,replace_dd$SSB/replace_dd$R)[60,]
cbind(replace_dd$SSB,replace_dd$R)[60,]

vpares$input$dat$maa
vpares$input$dat$waa

tmp = get.SPR(vpares)
plot(tmp$ysdata[-52,"SPR0"]~colSums(vpares$ssb)[-52])

plot(tmp$ysdata[-52,"SPR0"]~colSums(vpares$naa)[-52])

plot( (colSums(vpares$ssb)[-(51:52)]/1000),(colSums(vpares$ssb)[-(51:52)]/1000)/tmp$ysdata[-c(1,52),"SPR0"])

plot(replace_di$SSB,replace_di$R)
plot(replace_dd$SSB,replace_dd$R)






### calculating equilibrium ----


calc_eq = function(dd,SR,Fmulti=1) {
  
  if (isTRUE(dd)) {
    model_w0 = mod_w0
    model_wg = mod_w_growth
    model_mat = mod_mat
  } else {
    model_w0 = mod_w0_di
    model_wg = mod_w_growth_di
    model_mat = mod_mat_di
  }
  
  if(SR=="BH") resSR=resBH
  if(SR=="HS") resSR=resHS
  
  # Rc = seq(min(vpares$naa[1,]/1000)*0.01,0.5*max(as.numeric(vpares$naa[1,])/1000),length=200)
  
  M=0.4
  Fcurrent = rowMeans(vpares$faa[,as.character(2016:2020)])
  saa = Fcurrent/rev(Fcurrent)[1]
  
  faa = Fcurrent*Fmulti
  
  # l=50
  # 1:length(Rc) %>% map_dfr(., function(l) {

    # y <- 53.39
    obj_f2 = function(y,out=FALSE) {
      # summary(model_w0)
      weight = y # weight at age 0
      
      for (j in 1:A) {
        # if (j<A) {
        # summary(model_wg)
        preddata_wg = model_wg$data %>% 
          mutate(Age=j,Number_prev=n,Weight_prev=rev(weight)[1],Cohort_prev=nc[j],Cohort_plus=sum(nc[j:(j+1)]))
        weight = c(weight,as.numeric(predict(model_wg,newdata=preddata_wg[1,]))[1])
        # } else {
        #   n_A1 = r*exp(-(A-1)*M)
        #   n_plus = r*exp(-A*M)/(1-exp(-M))
        #   w_A1 = rev(weight)[1]
        #   # w_A1 = 682.14
        #   obj_f = function(x) {
        #     preddata_wg = model_wg$data %>% 
        #       mutate(Age=j,Number_prev=n,Weight_prev=(w_A1*n_A1+x*n_plus)/(n_A1+n_plus))
        #     tmp = as.numeric(predict(model_wg,newdata=preddata_wg[1,]))[1]
        #     return( (tmp-x)^2)
        #   }
        #   opt_w_plus = optimize(obj_f,interval=c(w_A1,2*w_A1))
        #   weight=c(weight,opt_w_plus$minimum)
        # }
      }
      # summary(model_mat)
      preddata_mat = model_mat$model %>% 
        mutate(Number_prev=n,Cohort_prev=nc[j],Cohort_plus=sum(nc[j:(j+1)]))
      maturity = c(0,sort(unique(as.numeric(predict(model_mat,newdata=preddata_mat)))),rep(1,3))
      ssb = sum(nc*weight*maturity)
      
      # summary(model_w0)
      preddata_w0=model_w0$data %>% mutate(Number_prev=n,Rec=nc[1],SSB=ssb)
      pred_w0 = as.numeric(predict(model_w0,newdata=preddata_w0[1,]))[1]
      
      if (out==FALSE) {
        return( (y-pred_w0)^2 )
      } else {
        res = list(naa=nc,waa=weight,maa=maturity,ssb=nc*weight*maturity)
        return( res )
      }
    }
    opt2 = optimize(obj_f2,interval=c(0.1,1000))
    # opt2$minimum
    # opt2$objective
    out = obj_f2(opt2$minimum,out=TRUE)
    out
    # out$ssb %>% sum
    
    refres = ref.F(Fcurrent=Fcurrent,Pope=TRUE,max.age=100,
                   maa=out$maa,waa=out$waa,M=rep(M,A+1),
                   waa.catch=out$waa,min.age=0,plot=FALSE)
    spr0 = refres$spr0
    
    # nc = sapply(0:(A-1), function(k) r*exp(-k*M))
    # nc = c(nc,n-sum(nc))
    # 
    ssb = sum(nc*weight*maturity)
    data.frame(SPR0=spr0,SSB=ssb,R=ssb/spr0)
  }
})



### weight and maturity regression ---- 


w0 <- glm(logW~1,data=dat)
w1 <- glm(logW~logN,data=dat)

AICc(w0,w1)
res_w = lrtest(w0,w1)


# m0 <- glm(logitM~1,data=dat)
# m1 <- glm(logitM~logN,data=dat)
# 
# summary(m1)
res_m = lrtest(m0,m1)

m0 = betareg(Maturity~1,data=dat,link="logit")
summary(m0)

m1 = betareg(Maturity~logN,data=dat,link="logit",type="BC")
summary(m1)

res_m = lrtest(m0,m1)

plot(dat$Maturity~predict(m1))

# N=100000
# tmp = runif(N,0,1)
# tmp_q = qnorm(tmp,0,1)
# hist(tmp_q)
# tmp_q %>% mean
# tmp_q %>% sd()

logit = function(p) log(p/(1-p))
inv_logit = function(x) 1/(1+exp(-x))

# # only latest five years
# w2 <- glm(logW~1,data=filter(dat,Year>max(Year)-5))
# w2$coefficients %>% exp()
# m2 <- glm(logitM~1,data=filter(dat,Year>max(Year)-5))
# m2$coefficients %>% inv_logit()

dat_wm = dat %>% 
  pivot_longer(values_to="Value",cols=c("Weight","Maturity"),
               names_to="Stat") %>%
  mutate(Stat=factor(Stat,levels=c("Weight","Maturity")))

qnorm(p=0.975)

pred_w = tibble(logN=seq(min(dat$logN),max(dat$logN),length=200))
tmp_w = allEffects(w1,xlevels=list(logN=pred_w$logN))
pred_w = pred_w %>%
  mutate(logN=seq(min(dat$logN),max(dat$logN),length=200),
         pred = tmp_w[[1]]$fit[,1],SE=tmp_w[[1]]$se) %>%
  mutate(upper=tmp_w[[1]]$upper[,1],lower=tmp_w[[1]]$lower[,1]) %>%
  mutate(N=exp(logN),Value=exp(pred),
         Upper=exp(upper),Lower=exp(lower),
         Stat = "Weight")

pred_m = tibble(logN=seq(min(dat$logN),max(dat$logN),length=200))
# tmp_m = predict(g1,newdata=pred_m,se=TRUE,type="link")
tmp_m = allEffects(m1,xlevels=list(logN=pred_m$logN))
# tmp_m[[1]]$fit
# tmp_m[[1]]$transformation$inverse(tmp_m[[1]]$fit)
# tmp_m[[1]]$se
# dat$Maturity
pred_m = tibble(logN=seq(min(dat$logN),max(dat$logN),length=200),
                pred = tmp_m[[1]]$fit[,1],
                SE=tmp_m[[1]]$se) %>%
  mutate(upper=tmp_m[[1]]$upper[,1],
         lower=tmp_m[[1]]$lower[,1]) %>%
  mutate(N=exp(logN),Value=inv_logit(pred),
         Upper=inv_logit(upper),Lower=inv_logit(lower),
         Stat = "Maturity")

pred_wm = bind_rows(pred_w,pred_m) %>%
  mutate(Stat=factor(Stat,levels=c("Weight","Maturity")))

base_size=12
point_size=1
path_size=1
g_wm = ggplot(data=dat_wm,aes(x=N,y=Value))+
  geom_ribbon(data=pred_wm,aes(ymin=Lower,ymax=Upper),alpha=0.4)+
  geom_path(data=pred_wm,size=path_size)+
  geom_point(size=point_size)+
  facet_wrap(vars(Stat),scales="free_y")+
  scale_x_log10()+
  # scale_y_log10()+
  theme_bw(base_size=base_size)+
  # ylim(0,NA)+
  xlab("Fish number (billion)")+ylab("")
g_wm

### Figure 1 ----

base_size=12
point_size=1.5
path_size=1

# weight figure 
dat_wm$Stat %>% unique()

g_w = ggplot(data=filter(dat_wm,Stat=="Weight"),aes(x=N,y=Value))+
  geom_ribbon(data=filter(pred_wm,Stat=="Weight"),aes(ymin=Lower,ymax=Upper),alpha=0.4)+
  geom_path(data=filter(pred_wm,Stat=="Weight"),size=path_size)+
  geom_point(size=point_size)+
  # facet_wrap(vars(Stat),scales="free_y")+
  scale_x_log10()+
  # scale_y_log10()+
  theme_bw(base_size=base_size)+
  # ylim(0,NA)+
  xlab("Fish number (billion)")+ylab("Body weight (g)")+
  ggtitle("(b)")
g_w


g_m = ggplot(data=filter(dat_wm,Stat=="Maturity"),aes(x=N,y=Value))+
  geom_ribbon(data=filter(pred_wm,Stat=="Maturity"),aes(ymin=Lower,ymax=Upper),alpha=0.4)+
  geom_path(data=filter(pred_wm,Stat=="Maturity"),size=path_size)+
  geom_point(size=point_size)+
  # facet_wrap(vars(Stat),scales="free_y")+
  scale_x_log10()+
  # scale_y_log10()+
  theme_bw(base_size=base_size)+
  # ylim(0,NA)+
  xlab("Fish number (billion)")+ylab("Maturity rate")+
  ggtitle("(c)")
g_m

(g_BH = ggplot(data=NULL,aes(x=SSB))+
    geom_ribbon(data=pred_BH,aes(ymax=upper,ymin=lower),alpha=0.4)+
    geom_path(data=pred_BH,aes(y=pred),size=path_size)+
    geom_point(data=SRdata,aes(y=R),size=point_size)+
    theme_bw(base_size=base_size)+
    xlab("SSB (thousand ton)")+ylab("Recruits (billion)") +
    ggtitle("(a)")
)


g_fig1 = gridExtra::grid.arrange(g_BH,g_w,g_m,nrow=1)


ggsave(g_fig1,filename="Fig_SR-weight-maturity.png",dpi=600,height=75,width=250,unit="mm")




### time-series figure ----

colnames(dat)

# dat2 = bind_rows(data.frame(Year=dat$Year,Value=dat$R,Stat="Fish number",Type="Age 0"),
#                  data.frame(Year=dat$Year,Value=dat$N,Stat="Fish number",Type="Older"),
#                  data.frame(Year=dat$Year,Value=dat$B,Stat="Biomass",Type="Total"),
#                  data.frame(Year=dat$Year,Value=dat$SSB,Stat="Biomass",Type="Spawner"),
#                  data.frame(Year=dat$Year,Value=dat$catch,Stat="Biomass",Type="Catch"),
#                  data.frame(Year=dat$Year,Value=dat$Weight,Stat="Weight",Type="Weight"),
#                  data.frame(Year=dat$Year,Value=dat$Maturity,Stat="Maturity",Type="Maturity"),
#                  data.frame(Year=dat$Year,Value=dat$F,Stat="F",Type="F")
#                  )


temp = frasyr::get.SPR(vpares)
ysdata = temp$ysdata[-nrow(temp$ysdata),]

dat = dat %>% bind_cols(ysdata)

plot(dat$SPR0/(dat$Weight*dat$Maturity))

plot(dat$SPR0~as.numeric(dat$Weight*dat$Maturity),log="xy")
plot(dat$SPR0/(dat$Weight*dat$Maturity)~dat$N,log="xy")

summary(glm(log(dat$SPR0/(dat$Weight*dat$Maturity))~dat$N))
summary(glm(log(dat$SPR0/(dat$Weight*dat$Maturity))~log(dat$N)))

c1 = glm(log(dat$SPR0)~1+offset(log(dat$Weight*dat$Maturity)))
library(effects)
plot(c1)

exp(mean(log(dat$SPR0/(dat$Weight*dat$Maturity))))

exp(-0.4)/(1-exp(-0.4))

plot((dat$SPR0/(dat$Weight*dat$Maturity))[-1]~dat$F[-50],log="y")

model = glm(log((dat$SPR0/(dat$Weight*dat$Maturity))[-1])~dat$F[-50])
summary(model)

model = glm(log((dat$SPR0/(dat$Weight*dat$Maturity))[-1])~dat$F[-50])
summary(model)

plot(dat$Weight,log(dat$Maturity/(1-dat$Maturity)),type="b")
plot(dat$Weight,dat$Maturity,type="b")


model=glm(log(dat$Maturity/(1-dat$Maturity))~dat$Weight)
summary(model)

dat2 = bind_rows(data.frame(Year=dat$Year,Value=dat$R,Stat="Fish number",Type="Age 0"),
                 data.frame(Year=dat$Year,Value=dat$N,Stat="Fish number",Type="Age 1+"),
                 data.frame(Year=dat$Year,Value=dat$B,Stat="Biomass",Type="Age 1+"),
                 data.frame(Year=dat$Year,Value=dat$SSB,Stat="Biomass",Type="Spawner"),
                 data.frame(Year=dat$Year,Value=dat$catch,Stat="Biomass",Type="Catch"),
                 data.frame(Year=dat$Year,Value=dat$Weight,Stat="Weight",Type="Weight"),
                 data.frame(Year=dat$Year,Value=dat$Maturity,Stat="Maturity",Type="Maturity"),
                 data.frame(Year=dat$Year,Value=dat$F,Stat="F",Type="F")
)

head(dat2)

# SPR0


plot(spr0~dat$N,log="")
plot(spr0~dat$N,log="x")
plot(spr0~dat$N,log="xy")

spr0_res= glm(log(spr0)~dat$logN)
summary(spr0_res)

spr0_res2= glm(log(spr0)~dat$N)
summary(spr0_res2)

AICc(spr0_res,spr0_res2)

spr0_2 = 1:nrow(dat) %>% map_dbl(function(i) {
  temp = frasyr::calc_steepness(SR="BH",rec_pars=resBH$pars,M=0.4,
                             waa=c(0,dat$Weight[i]), maa=c(0,dat$Maturity[i]),
                       plus_group = FALSE)
  temp[1,"SPR0"]
})

plot(spr0_2~dat$N,log="x")
plot(spr0_2~dat$N,log="xy")


matplot(cbind(spr0,spr0_2))

plot(spr0)

path_size=0.8
base_size=10

# theme()

g1 = ggplot(data=filter(dat2,Stat=="Fish number"),aes(x=Year,y=Value))+
  geom_path(aes(colour=Type,linetype=Type),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="",guide=guide_legend(nrow=1))+
  scale_linetype_discrete(name="",guide=guide_legend(nrow=1))+
  ylim(0,NA)+labs(colour=NULL,linetype=NULL)+
  ylab("Number (billion)")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(-2, "pt"))
g1

g2 = ggplot(data=filter(dat2,Stat=="Biomass") %>% 
              mutate(Type=factor(Type,levels=c("Age 1+","Spawner","Catch"))),aes(x=Year,y=Value))+
  geom_path(aes(colour=Type,linetype=Type),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="",guide=guide_legend(nrow=1))+
  scale_linetype_discrete(name="",guide=guide_legend(nrow=1))+
  ylim(0,NA)+labs(colour=NULL,linetype=NULL)+
  ylab("Biomass (1000 ton)")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(-2, "pt"))
g2

g3 = ggplot(data=filter(dat2,Stat=="Weight"),aes(x=Year,y=Value))+
  # geom_path(aes(colour=Type,linetype=Type),size=path_size)+
  geom_path(size=path_size)+
  theme_bw(base_size=base_size)+
  # scale_colour_brewer(palette="Set1",name="")+
  # scale_linetype_discrete(name="")+
  ylim(0,NA)+
  ylab("Body weight (g)")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1))
g3

g4 = ggplot(data=filter(dat2,Stat=="Maturity"),aes(x=Year,y=Value))+
  # geom_path(aes(colour=Type,linetype=Type),size=path_size)+
  geom_path(size=path_size)+
  theme_bw(base_size=base_size)+
  # scale_colour_brewer(palette="Set1",name="")+
  # scale_linetype_discrete(name="")+
  ylim(0,NA)+
  ylab("Maturation rate")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1))
g4

g5 = ggplot(data=filter(dat2,Stat=="F"),aes(x=Year,y=Value))+
  # geom_path(aes(colour=Type,linetype=Type),size=path_size)+
  geom_path(size=path_size)+
  theme_bw(base_size=base_size)+
  # scale_colour_brewer(palette="Set1",name="")+
  # scale_linetype_discrete(name="")+
  ylim(0,NA)+
  ylab("F")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1))
g5

g_time = gridExtra::grid.arrange(g1,g2,g3,g4,g5,nrow=2)

ggsave(g_time,filename="graph_timeseries.png",dpi=600,
       unit="mm",height=120,width=240)


# pred_dat$SSB %>% range
pred_dat = data.frame(logN = seq(min(dat$logN)-1,max(dat$logN)+0.1,length=1000)) %>%
  mutate(N = exp(logN),
         Weight=exp(predict(w1,newdata=.)),
         Maturity=inv_logit(predict(m1,newdata=.))) %>%
  mutate(B = Weight*N) %>%
  mutate(SSB = Maturity*B,Density="Dependent")

pred_dat2 = data.frame(logN = seq(min(dat$logN)-1,max(dat$logN)+0.1,length=1000)) %>%
  mutate(N = exp(logN),
         Weight=exp(predict(w0,newdata=.)),
         Maturity=inv_logit(predict(m0,newdata=.))) %>%
  mutate(B = Weight*N) %>%
  mutate(SSB = Maturity*B,Density="Independent")

# pred_dat3 = data.frame(logN = seq(min(dat$logN)-1,max(dat$logN)+0.1,length=1000)) %>%
#   mutate(N = exp(logN),
#          Weight=exp(predict(w2,newdata=.)),
#          Maturity=inv_logit(predict(m2,newdata=.))) %>%
#   mutate(B = Weight*N) %>%
#   mutate(SSB = Maturity*B,Type="DI_latest5")

# pred_dat = bind_rows(pred_dat,pred_dat2,pred_dat3)

pred_dat = bind_rows(pred_dat,pred_dat2)

rec_pars_BH = as.list(resBH$pars)
rec_pars_HS = as.list(resHS$pars)


# i=1
# calc_steepness(SR="HS",rec_pars=rec_pars_HS,plus_group = TRUE,
#                waa = c(0,pred_dat$Weight[i]),
#                maa = c(0,pred_dat$Maturity[i]),
#                M = rep(0.4,2))
# 
# calc_steepness(SR="BH",rec_pars=rec_pars_BH,plus_group = TRUE,
#                waa = c(0,pred_dat$Weight[i]),
#                maa = c(0,pred_dat$Maturity[i]),
#                M = rep(0.4,2))

pred_dat2 = 1:nrow(pred_dat) %>% map_dfr(., function(i) {
  bind_cols(pred_dat[i,],calc_steepness(SR="BH",rec_pars=rec_pars_BH,plus_group = TRUE,
                                        waa = c(0,pred_dat$Weight[i]),
                                        maa = c(0,pred_dat$Maturity[i]),
                                        M = rep(0.4,2))
  )
})

# pred_dat2$Type %>% unique()

# pred_dat2 %>% filter(Type=="DI_latest5")

pred_dat2 = pred_dat2 %>% 
  mutate(R = SSB/SPR0) %>%
  mutate(R_SR = SSB %>% map_dbl(., function(x) SRF_BH(x,a=rec_pars_BH$a,b=rec_pars_BH$b))) %>%
  mutate(SSB_excess = R_SR*SPR0 ) %>%
  mutate(SPS = SSB_excess/SSB)


# plot(R ~ SSB,data=pred_dat2)

# data_line =select(pred_dat2,R,SSB) %>% mutate(Type="Replacement") %>%
#   # full_join(resBH$pred %>% mutate(Type="Beverton-Holt")) %>%
#   # full_join(resHS$pred %>% mutate(Type="Hockey-Stick")) %>%
#   mutate(Type = factor(Type,levels=c("Replacement","Beverton-Holt","Hockey-Stick")))

g1_BH = ggplot(data=NULL,aes(x=SSB,y=R))+
  geom_path(data=resBH$pred,size=path_size,colour="black",linetype="longdash")+
  geom_point(data=as.data.frame(SRdata))+
  geom_path(data=filter(pred_dat2,SSB<max(resBH$pred$SSB)),aes(colour=Density,linetype=Density),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="Density")+
  scale_linetype_discrete(name="Density")+
  xlab("Spawning stock biomass (1000 ton)")+
  ylab("Number of recruits (billion)")+
  theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
        legend.margin = margin(0, 0, 0, 0))

g1_BH

ggsave(g1_BH,filename="SR_replacement_BH.png",dpi=600,units="mm",
       height=80,width=120)


# g1_HS = ggplot(data=NULL,aes(x=SSB,y=R))+
#   geom_path(data=resHS$pred,size=path_size,colour="black",linetype="longdash")+
#   geom_point(data=as.data.frame(SRdata))+
#   geom_path(data=filter(pred_dat2,SSB<max(resBH$pred$SSB)),aes(colour=Density,linetype=Density),size=path_size)+
#   theme_bw(base_size=base_size)+
#   scale_colour_brewer(palette="Set1",name="Density")+
#   scale_linetype_discrete(name="Density")+
#   xlab("Spawning stock biomass (1000 ton)")+
#   ylab("Number of recruits (billion)")+
#   theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
#         legend.margin = margin(0, 0, 0, 0))
# 
# g1_HS
# 
# ggsave(g1_HS,filename="SR_replacement_HS.png",dpi=600,units="mm",
#        height=80,width=120)

g1_SPS = ggplot(data=NULL,aes(x=SSB,y=SPS))+
  geom_path(data=filter(pred_dat2,SSB<max(resBH$pred$SSB)),aes(colour=Density,linetype=Density),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="Density")+
  scale_linetype_discrete(name="Density")+
  xlab("Spawning stock biomass (1000 ton)")+
  ylab("Spawners per spwaner")+
  theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
        legend.margin = margin(0, 0, 0, 0))

g1_SPS

ggsave(g1_SPS,filename="SPS_BH.png",dpi=600,units="mm",
       height=80,width=120)

write.csv(pred_dat2, file="SPR_SPS.csv",row.names=FALSE)

pred_dat[100,]

# F = 0.0001
# faa = c(0,f)
# n = 10
obj_fun = function(n,F,model_w=w1,model_m=m1,resSR=resBH,out=TRUE) {
  f = F
  tmp = data.frame(logN=log(n))
  m_pred = as.numeric(inv_logit(predict(model_m,newdata=tmp)[1]))
  w_pred = as.numeric(exp(predict(model_w,newdata=tmp)[1]))
  refres = ref.F(Fcurrent=c(0,f),Pope=TRUE,max.age=100,
                 maa=c(0,m_pred),waa=c(0,w_pred),M=c(0.4,0.4),
                 waa.catch=c(0,w_pred),min.age=0,plot=FALSE)
  SPR = refres$currentSPR$SPR
  SSB = as.numeric(n*m_pred*w_pred)
  a = resSR$pars$a
  b = resSR$pars$b
  SR = resSR$input$SR

  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  
  R_pred1 = SRF(x=SSB,a=a,b=b)
  R_pred2 = SSB/SPR
  if (out==TRUE) {
    (R_pred1-R_pred2)^2
  } else {
    Yield = n*w_pred*exp(-0.4/2)*(1-exp(-f))
    tibble(N=n,Weight=w_pred,Maturity=m_pred,B=n*w_pred,SSB=SSB,F=f,SPR=SPR,R=R_pred1,Yield=Yield,U=Yield/B)
  }
}

range(dat$N)
x = 0.00001
# f = seq(0.0001,100,length=1000)
f = c(exp(seq(log(0.001),log(10),length=901)))
# x = f[2]
dd_f = f %>% map_dfr(., function(x) {

  # tmp_n = exp(seq(log(0.1),log(100),length=31))
  # tmp_obj = tmp_n %>% map_dbl(., function(i) obj_fun(n=i,F=x,
  #                   model_w=w1,model_m=m1,resSR=resBH,out=TRUE))
  # int = tmp_n[c(max(1,which.min(tmp_obj)-1),min(length(tmp_obj),which.min(tmp_obj)+1))]
  int = c(0.1,1000)
  opt = optimize(obj_fun, interval = int,F=x,
                 model_w=w1,model_m=m1,resSR=resBH,out=TRUE)

  tbl = obj_fun(opt$minimum,F=x,model_w=w1,model_m=m1,resSR=resBH,out=FALSE)
  tbl
})

# x = 0.3

opt_msy = function(x,out=TRUE,model_w=w1,model_m=m1,resSR=resBH) {
  int = c(0.1,1000)
  opt = optimize(obj_fun, interval = int,F=x,
                 model_w=model_w,model_m=model_m,resSR=resSR,out=TRUE)
  tbl = obj_fun(opt$minimum,F=x,model_w=model_w,model_m=model_m,resSR=resSR,out=FALSE)
  if (out==TRUE) {
    -pull(tbl,Yield)[1]
  } else {
    tbl
  }
}

msy_dd_opt = optimize(opt_msy,interval=c(0.1,2),model_w=w1,model_m=m1,resSR=resBH)
msy_dd = opt_msy(msy_dd_opt$minimum,model_w=w1,model_m=m1,resSR=resBH,out=FALSE)

opt_msy(10,model_w=w0,model_m=m0,resSR=resBH)

msy_di_opt = optimize(opt_msy,interval=c(0.1,1),model_w=w0,model_m=m0,resSR=resBH)
msy_di = opt_msy(msy_di_opt$minimum,model_w=w0,model_m=m0,resSR=resBH,out=FALSE)

msy_data = full_join(
  msy_dd %>% mutate(Density="Dependent"),
  msy_di %>% mutate(Density="Independent")
)

write.csv(msy_data,file="msy_data.csv")

plot(Yield~SSB,data=dd_f,type="l")
plot(Yield/SSB~SSB,data=dd_f,type="l",log="")

plot(Yield~F,data=dd_f,type="l")
dd_f$Yield
plot(catch/SSB~SSB,data=dat,type="p",log="")

# c(tbl$R,tbl$SSB/tbl$SPR)

f2 = c(exp(seq(log(0.001),log(10),length=901)))
# x=0.64
# f2 = c(seq(0.62,0.64,length=11))

di_f = f2 %>% map_dfr(., function(x) {
  
  # tmp_n = exp(seq(log(0.1),log(100),length=31))
  # tmp_obj = tmp_n %>% map_dbl(., function(i) obj_fun(n=i,F=x,
  #                   model_w=w1,model_m=m1,resSR=resBH,out=TRUE))
  # int = tmp_n[c(max(1,which.min(tmp_obj)-1),min(length(tmp_obj),which.min(tmp_obj)+1))]
  int = c(0.1,1000)
  opt = optimize(obj_fun, interval = int,F=x,
                 model_w=w0,model_m=m0,resSR=resBH,out=TRUE)
  
  tbl = obj_fun(opt$minimum,F=x,model_w=w0,model_m=m0,resSR=resBH,out=FALSE)
  tbl
})

di_f

plot(Yield~SSB,data=di_f,type="l")
plot(Yield/SSB~SSB,data=di_f,type="l",log="")

plot(Yield~F,data=dd_f,type="l",log="")
plot(Yield~F,data=di_f,type="l",log="")

plot(SSB~F,data=di_f,type="l",log="")

plot(SSB~F,data=dd_f,type="l",log="")

plot(catch/SSB~SSB,data=dat,type="p",log="")

d_f = full_join(
  dd_f %>% mutate(Density="Dependent"),
  di_f %>% mutate(Density="Independent")) %>%
  mutate(YPS = Yield/SSB)

di_f

write.csv(d_f, file="Sustainable_yield.csv",row.names=FALSE)


g_yc = ggplot(filter(d_f,SSB<=max(dd_f$SSB) & SSB>=min(dd_f$SSB)),aes(x=SSB,y=Yield))+
  geom_path(size=path_size,aes(colour=Density,linetype=Density))+
  # geom_point(data=dat,aes(x=SSB,y=catch))+
  theme_bw(base_size=base_size)+
  scale_color_brewer(palette="Set1")+
  ylab("Sustainable yield (1000 ton)")+
  xlab("Spawning stock biomass (1000 ton)")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt"))

g_yc

g_surplus = ggplot(filter(d_f,SSB<=max(dd_f$SSB) & SSB>=min(dd_f$SSB)),aes(x=SSB,y=YPS))+
  geom_path(size=path_size,aes(colour=Density,linetype=Density))+
  # geom_point(data=dat,aes(x=SSB,y=catch))+
  theme_bw(base_size=base_size)+
  scale_color_brewer(palette="Set1")+
  ylab("Sustainable yield relative to SSB")+
  xlab("Spawning stock biomass (1000 ton)")+
  theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt"))

g_surplus

head(d_f)

d_YS = d_f %>% 
  mutate(MSY = if_else(Density=="Dependent",msy_dd$Yield,msy_di$Yield),
         SSB0 = if_else(Density=="Dependent",max(dd_f$SSB),max(di_f$SSB))
  ) %>%
  mutate("SY/MSY" = Yield/MSY,"SSB/SSB0"=SSB/SSB0) %>%
  pivot_longer(cols=c("SY/MSY","SSB/SSB0"),names_to="Reference",values_to="Value") %>%
  mutate(Reference = factor(Reference,levels=c("SY/MSY","SSB/SSB0")))

# d_YS$`SY/MSY` %>% range
# d_YS$`SSB/SSB0` %>% range

g_F = ggplot(filter(d_YS,F<1.5),aes(x=F,y=Value))+
  geom_path(size=path_size,aes(colour=Density,linetype=Reference))+
  # geom_point(data=dat,aes(x=SSB,y=catch))+
  theme_bw(base_size=base_size)+
  scale_color_brewer(palette="Set1")+
  scale_linetype_discrete(name="Type")+
  ylab("Relative value")+
  xlab("Fishing mortality coefficient")+
  theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt"),
        legend.key.size=unit(0.4,"cm"))+
  guides(colour=guide_legend(nrow=1,order=1,title.position="top"),
         linetype=guide_legend(nrow=1,order=2,title.position="top"))+
  ylim(0,1.25)

g_F



g1_BH = ggplot(data=NULL,aes(x=SSB,y=R))+
  geom_path(data=resBH$pred,size=path_size,colour="black",linetype="longdash")+
  geom_point(data=as.data.frame(SRdata))+
  geom_path(data=filter(pred_dat2,SSB<max(resBH$pred$SSB)),aes(colour=Density,linetype=Density),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="Density")+
  # scale_linetype_manual(name="")+
  scale_linetype_discrete(name="Density")+
  xlab("Spawning stock biomass (1000 ton)")+
  ylab("Number of recruits (billion)")+
  theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt"))
# +
#   geom_vline(xintercept=msy_dd$SSB,colour="orange",linetype="dotted",size=path_size)

g1_BH

ggsave(g1_BH,filename="SR_replacement_BH_noMSY.png",dpi=600,units="mm",
       height=80,width=120)

# g1_BH <- g1_BH+ggtitle("(a)")

# g_2 = gridExtra::grid.arrange(g1_BH+ggtitle("(a)"),
#                               g_yc+ggtitle("(b)"),
#                               g_surplus+ggtitle("(c)"),
#                               g_F+ggtitle("(d)"),
#                               ncol=2)

g1_SPS = ggplot(data=NULL,aes(x=SSB,y=SPS))+
  geom_path(data=filter(pred_dat2,SSB<max(resBH$pred$SSB)),aes(colour=Density,linetype=Density),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="Density")+
  scale_linetype_discrete(name="Density")+
  xlab("Spawning stock biomass (1000 ton)")+
  ylab("Spawners per spwaner")+
  theme(legend.position=c(0.99,0.99),legend.justification=c(1,1),
        legend.margin = margin(0, 0, 0, 0))+
  # geom_hline(yintercept=1,linetype="dotted",colour="darkgray",size=path_size)+
  ylim(0,NA)

g1_SPS


g_yc2 = g_yc +
    geom_vline(xintercept=msy_dd$SSB,colour="darkorange",linetype="dotted",size=path_size)
g_yc2

  
ggsave(g_yc2,filename="YieldCurve_BH.png",dpi=600,units="mm",
       height=80,width=120)


ggsave(g_F,filename="F_BH.png",dpi=600,units="mm",
       height=80,width=120)



g_2 = gridExtra::grid.arrange(g1_BH+ggtitle("(a)"),
                              g1_SPS+ggtitle("(b)"),
                              g_yc+ggtitle("(c)"),
                              # g_surplus+ggtitle("(c)"),
                              g_F+ggtitle("(d)"),
                              ncol=2)


ggsave(g_2,filename="Population-level_density-dependence.png",dpi=600,units="mm",
       height=125,width=200)



### time series figures with MSY ref ----

head(dat2)

dat3 = dat2 %>%
  mutate(Category = "Annual value")


dat3$Stat %>% unique()
dat3$Type %>% unique()

msy_longer = msy_data %>% filter(Density=="Dependent") %>%
  pivot_longer(names_to = "stat", values_to="Value",cols=-Density) %>%
  mutate(Stat=case_when(stat=="N" |stat=="R" ~ "Fish number",
                        stat=="B" |stat=="SSB" | stat=="Yield" ~ "Biomass",
                        TRUE ~ stat)) %>%
  mutate(Type=case_when(stat=="R" ~ "Age 0",
                        stat=="N" ~ "Older",
                        stat=="B" ~ "Total",
                        stat=="SSB" ~ "Spawner",
                        stat=="Yield" ~ "Catch",
                        TRUE ~ stat)) %>%
  mutate(Category = "MSY level")

dat3 = full_join(dat3,msy_longer)


g1 = ggplot(data=filter(dat2,Stat=="Fish number"),aes(x=Year,y=Value))+
  geom_path(aes(colour=Type),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="",guide=guide_legend(nrow=1))+
  # scale_linetype_discrete(name="",guide=guide_legend(nrow=1))+
  ylim(0,NA)+labs(colour=NULL,linetype=NULL)+
  ylab("Number (billion)")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(-2, "pt"))+
  geom_hline(data=filter(msy_longer,Stat=="Fish number"),aes(yintercept=Value,colour=Type),
             size=0.5,linetype="dotted")
g1

g2 = ggplot(data=filter(dat2,Stat=="Biomass") %>% mutate(Type=factor(Type,levels=c("Total","Spawner","Catch"))),aes(x=Year,y=Value))+
  geom_path(aes(colour=Type),size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="",guide=guide_legend(nrow=1))+
  # scale_linetype_discrete(name="",guide=guide_legend(nrow=1))+
  ylim(0,NA)+labs(colour=NULL,linetype=NULL)+
  ylab("Biomass (1000 ton)")+
  theme(legend.position=c(0.01,0.99),legend.justification=c(0,1),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "pt"),
        legend.spacing.y = unit(-2, "pt"))+
  geom_hline(data=filter(msy_longer,Stat=="Biomass"),aes(yintercept=Value,colour=Type),
             size=0.5,linetype="dotted")

g2

g3 = ggplot(data=filter(dat2,Stat=="Weight"),aes(x=Year,y=Value))+
  # geom_path(aes(colour=Type,linetype=Type),size=path_size)+
  geom_path(size=path_size,aes(colour=Type))+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="")+
  # scale_linetype_discrete(name="")+
  ylim(0,NA)+
  ylab("Body weight (g)")+
  theme(legend.position="",legend.justification=c(0,1))+
  geom_hline(data=filter(msy_longer,Stat=="Weight"),aes(yintercept=Value,colour=Type),
             size=0.5,linetype="dotted")

g3

g4 = ggplot(data=filter(dat2,Stat=="Maturity"),aes(x=Year,y=Value))+
  geom_path(aes(colour=Type),size=path_size)+
  # geom_path(size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="")+
  # scale_linetype_discrete(name="")+
  ylim(0,NA)+
  ylab("Maturation rate")+
  theme(legend.position="",legend.justification=c(0,1))+
  geom_hline(data=filter(msy_longer,Stat=="Maturity"),aes(yintercept=Value,colour=Type),
             size=0.5,linetype="dotted")

g4

g5 = ggplot(data=filter(dat2,Stat=="F"),aes(x=Year,y=Value))+
  geom_path(aes(colour=Type),size=path_size)+
  # geom_path(size=path_size)+
  theme_bw(base_size=base_size)+
  scale_colour_brewer(palette="Set1",name="")+
  # scale_linetype_discrete(name="")+
  ylim(0,NA)+
  ylab("F")+
  theme(legend.position="",legend.justification=c(0,1))+
  geom_hline(data=filter(msy_longer,Stat=="F"),aes(yintercept=Value,colour=Type),
             size=0.5,linetype="dotted")

g5

g_time = gridExtra::grid.arrange(g1+ggtitle("(a)"),
                                 g2+ggtitle("(b)"),
                                 g3+ggtitle("(c)"),
                                 g4+ggtitle("(d)"),
                                 g5+ggtitle("(e)"),nrow=2)

ggsave(g_time,filename="graph_timeseries_MSY.png",dpi=600,
       unit="mm",height=120,width=240)


### FishLife ----

library(FishLife)

( Predictions = Plot_taxa(Search_species(Genus="Scomber",Species="japonicus")$match_taxonomy) )
Predictions[[1]]$Mean_pred["ln_MASPS"] %>% exp

# あんまり合わない





# pred_dat = data.frame(N = seq(exp(min(dat$logN)-1),exp(max(dat$logN)+1),length=1000)) %>%
#   mutate(logN = log(N)) %>%
#   mutate(Weight=exp(predict(m3,newdata=.)),
#          Maturity=inv_logit(predict(l3,newdata=.))) %>%
#   mutate(B = Weight*N) %>%
#   mutate(SSB = Maturity*B)

plot(Weight~N,data=pred_dat)
plot(Maturity~N,data=pred_dat)

plot(Weight~SSB,data=pred_dat)
plot(Maturity~SSB,data=pred_dat)

plot(pred_dat$SSB~pred_dat$N)

plot(dat$SSB~dat$N)


range(dat$N)
range(dat$SSB)
as.numeric(logLik(m1))

head(dat)

plot(dat)

summary(glm(logW ~ N,data=dat))
summary(glm(logW ~ log(N),data=dat))

model1 = glm(logW ~ log(N),data=dat)

plot(model1)

summary(glm(logitM ~ N,data=dat))
summary(glm(logitM ~ log(N),data=dat))

model2 = glm(logitM ~ log(N),data=dat)
plot(model2)

summary(glm(logitM ~ B,data=dat))
summary(glm(logitM ~ log(B),data=dat))

logit = function(p) log(p/(1-p))
inv_logit = function(y) 1/(1+exp(-y))

mu = -1
sigma = 1
obj = function(p,mu,sigma) {
    tmp = integrate(function(x) inv_logit(x)*dnorm(x,mean=logit(p),sd=sigma),lower=-Inf,upper=+Inf)
    (tmp$value - inv_logit(mu))^2
}


opt = optimize(obj,c(0.000001,0.999999),mu=mu,sigma=sigma)


N = 10000
y = rnorm(N,-1,sd=1)
inv_logit(-1)
p = inv_logit(y)

library(nlme)

modelW0 = glm(log(Weight)~log(N),data=dat)
modelW1 = glm(log(Weight)~N,data=dat)
AICc(modelW0,modelW1)


modelMat0 = glm(logit(Maturity)~log(N),data=dat)
summary(modelMat0)

modelMat1 = glm(logit(Maturity)~N,data=dat)

# AICc(modelMat0,modelMat1)

plot(modelMat0$residuals)
stats::ar(modelMat0$residuals,order.max = 1)

arima_res = arima(dat[,"logW"],order=c(0,0,0),
                  method="ML",xreg=dat[,"logN"])
arima_res$coef
AICc(arima_res)

arima_res = arima(dat[,"logW"],order=c(1,0,0),
                  method="ML",xreg=dat[,"logN"])
arima_res$coef
AICc(arima_res)



AICc(m)
dat$t

# install.packages("brms")
library(brms)
d <- simCor1(phi=0.8,sdgrp=2,sdres=1,seed=101)



modelW0_N = glm(log(Weight)~N,data=dat)
modelW_RE_N = lme(log(Weight)~N,random=~1|YEAR,data=dat)
modelW_AR_N = lme(log(Weight)~N,random=~1|YEAR,data=dat,correlation=corAR1())


summary(modelW_AR)
AIC(modelW0,modelW_RE,modelW_AR,modelW0_N,modelW_RE_N,modelW_AR_N)

logLik(modelW0)
logLik(modelW_AR)

plot(modelW0$residuals,type="l")

summary(modelW_AR)

modelW0

library(MuMIn)

modelW_RE2 = glmmTMB(log(Weight)~log(N)+(1|YEAR),data=dat)
summary(modelW_RE2)

AICtab(modelW_RE2)
AICc(modelW_RE2)

modelW = lme(log(Weight)~log(N),random=~1|YEAR,data=dat,correlation=corAR1())
summary(modelW)
AIC(modelW)


?simCor1
?glmmTMB

log(1)


summary(p)
