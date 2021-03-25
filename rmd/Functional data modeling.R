library(readr)
library(fdapace)
library(tidyverse)
library(plyr)
library(dplyr)
library(fda)
library(lubridate)
library(ggplot2)
library(fdapace)
library(splines)
library(Matrix)
library(fds)
library(rainbow)
library(MASS)
library(pcaPP)
library(RCurl)
library(fda.usc)
library(mboost)
library(parallel)
library(stabs)
library(FDboost)

##### two years data - data__all ----------------------------------------------

# outlier -----------------------------------------------------------------
hour <- read.csv("hour.csv")
# filter
data_all0_all <- hour %>% mutate(date=ymd(dteday)) %>%  dplyr::select(hr,cnt,date)  
data_all <- data_all0_all %>% pivot_wider(
  names_from = date,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.data.frame()
# set up
data_all1 <- t(data_all) [-1,]%>% as.data.frame()
data_all_fdata_all <- fdata(data_all1)
plot(optim.np(data_all_fdata_all)$fdata_all.est)
# detecting outliers
outliers.depth.pond(data_all_fdata_all,nb=200) # outliers "2012-07-04" "2012-09-08" # dep.out 9.707107 9.336198
# plotting
data0 <- hour %>% mutate(date=ymd(dteday)) %>%  dplyr::select(hr,cnt,date)  
data <- data0 %>% pivot_wider(
  names_from = date,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.data.frame()
# set up
data1 <- t(data) [-1,]%>% as.data.frame()
data_fdata <- fdata(data1)
plot(optim.np(data_fdata)$fdata.est)
# detecting outliers
outliers.depth.pond(data_fdata,nb=200) # 18,26,27
# plot outliers
plot(data_fdata,col="grey",main="Outlier detection for data in Jan 2011",
     xlab = "time [hour]", ylab = "count" )

out1 <- data0 %>% filter(date== "2012-07-04"| date=="2012-09-08")  
out2 <- out1 %>% pivot_wider(
  names_from = date,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.data.frame()
out <- t(out2) [-1,]%>% as.data.frame()
data_fdata_out <- fdata(out)
plot(data_fdata,col="grey",main="Outlier detection for data in year 2011 and 2012",
     xlab = "time [hour]", ylab = "count" )
lines(data_fdata_out ,lwd=2,lty=2:4,col=2:4)

# set up dataset ----------------------------------------------------------
# filter
data_all0_all <- hour %>% mutate(date=ymd(dteday))%>% 
  filter(date %in% c("2012-07-04","2012-09-08")== FALSE)
# count
count <- data_all0_all %>%  dplyr::select(hr,cnt,date) %>%
  mutate(hr=paste("hour",hr))  %>% pivot_wider(
  names_from = date,
  names_prefix = "date_",
  values_from = cnt,
  values_fill = list(cnt=0)
) %>% as.data.frame() %>% dplyr::select(-hr) 
# temp
temp <- data_all0_all %>%  dplyr::select(hr,temp,date) %>% 
  mutate(hr=paste("hour",hr))   %>% pivot_wider(
  names_from = date,
  names_prefix = "date_",
  values_from = temp,
  values_fill = list(temp=0)
) %>% as.data.frame() %>% dplyr::select(-hr) 
# work
work <- data_all0_all %>%  dplyr::select(hr,workingday,date) %>%
  mutate(hr=paste("hour",hr))%>% pivot_wider(
  names_from = date,
  names_prefix = "date_", 
  values_from = workingday,
  values_fill = list(workingday=0)
)  %>% filter(hr=="hour 1") %>% dplyr::select(-hr) %>% as.matrix()
work <- t(work) %>% as.data.frame() %>% mutate(work=as.factor(V1)) %>%  mutate(work = work 
                                        %>%fct_recode(
                                        "neither weekend nor holidate" = "1",
                                        "weekend nor holidate" = "0"
                                      )) %>% dplyr::select(work)  %>% as.matrix()
# wind
wind <- data_all0_all %>%  dplyr::select(hr,windspeed,date) %>% 
  mutate(hr=paste("hour",hr))  %>% pivot_wider(
  names_from = date,
  names_prefix = "date_",
  values_from = windspeed,
  values_fill = list(windspeed=0)
) %>% as.data.frame() %>% dplyr::select(-hr) 
# week
week <- data_all0_all %>%  dplyr::select(hr,weekday,date)  %>%
  mutate(hr=paste("hour",hr))  %>% pivot_wider(
  names_from = date,
  names_prefix = "date_",
  values_from = weekday,
  values_fill = list(week=0)
) %>% filter(hr=="hour 1") %>% dplyr::select(-hr) %>% as.matrix()
week <- t(week) %>% 
  as.data.frame() %>% mutate(week=as.factor(V1)) %>% 
  mutate(week = week %>% fct_recode( "Sundate" = "0",
                                     "Mondate" = "1",
                                     "Tuesdate" = "2",
                                     "Wednesdate" = "3",
                                     "Thursdate" = "4",
                                     "Fridate" = "5",
                                     "Saturdate"= "6")) %>% dplyr::select(week) %>% 
  mutate(week=week %>% fct_relevel("Saturdate", "Sundate"))  %>% drop_na() %>% as.matrix()
# date
date <- data_all0_all %>% dplyr::select(hr,weekday,date)%>% 
  mutate(date=paste("date",date)) %>% pivot_wider(
  names_from = hr,
  names_prefix = "hr_",
  values_from = weekday,
  values_fill = list(week=0) 
)  %>% dplyr::select(date) %>% mutate(date=as.factor(date)) %>% as.matrix()
# weather
weather <- data_all0_all %>%  dplyr::select(hr,weathersit,date)  %>% 
  mutate(hr=paste("hour",hr))  %>% pivot_wider(
  names_from = date,
  names_prefix = "date_",
  values_from = weathersit,
  values_fill = list(weathersit=1)
) %>% dplyr::select(-hr)  %>% as.matrix()
# hum
humidity <- data_all0_all %>%  dplyr::select(hr,hum,date)  %>% 
  mutate(hr=paste("hour",hr))  %>% pivot_wider(
  names_from = date,
  names_prefix = "date_",
  values_from = hum,
  values_fill = list(hum=0)
) %>% dplyr::select(-hr)
# season
season <- data_all0_all %>%  dplyr::select(hr,season,date)  %>% 
  mutate(hr=paste("hour",hr))  %>% pivot_wider(
  names_from = date,
  names_prefix = "date_", 
  values_from = season,
  values_fill = list(season=3)
)  %>% filter(hr=="hour 1") %>% dplyr::select(-hr) %>% as.matrix()
season <- t(season) %>% as.data.frame() %>% 
  mutate(season=as.factor(V1)) %>%  mutate(season = season %>% 
                                             fct_recode("winter" = "1",
                                                        "spring" = "2",
                                                        "summer" = "3",
                                                         "fall" = "4")) %>% 
  dplyr::select(season) %>% as.matrix()

# data__all
data__all <- list(count = t(count), temp = t(temp),wind=t(wind), 
                  work=as.factor(work),week=as.factor(week), 
                  date=as.factor(date),
                  humidity=t(humidity),season=as.factor(season))

# define function indices
data__all$hour.t <- 0:23
data__all$hour.s <- 0:23

# center temperature curves:
# data__all$temp <- sweep(data__all$temp, 2, colMeans(data__all$temp))
data__all$humidity <- scale(data__all$humidity,scale = F)
data__all$temp <- scale(data__all$temp,scale = F)
data__all$weather <- scale(data__all$weather,scale = F)
data__all$wind <- scale(data__all$wind,scale = F)

# k min -------------------------------------------------------------------

k_min10 <- FDboost(count ~ 1, timeformula = ~ bbs(hour.t, knots = 11, cyclic = TRUE, df=3), offset_control = o_control(k_min = 10), data=data__all)
k_min15 <- FDboost(count ~ 1, timeformula = ~ bbs(hour.t, knots = 11, cyclic = TRUE, df=3), offset_control = o_control(k_min = 15), data=data__all)

par(mfrow=c(2,2))
plot(k_min10,ask=F,main="offset k_min=10")
plot(k_min15,ask=F,main="offset k_min=15") # choose K_min=15

# final model ----------------------------------------------------------------
final <- FDboost(count ~ 1 
                 + bsignal(temp, hour.s, knots = 10, df = 4)
                 + bsignal(humidity, hour.s, knots = 10,  df = 4)
                 + bols(work,df=1), 
                 timeformula = ~ bbs(hour.t, knots = 10,  df=4),  
                 offset_control = o_control(k_min = 15), 
                 control=boost_control(mstop = 100,nu=0.001),
                 family=Poisson(),
                 data=data__all)
stabsel(final,q=4, PFER = 1)
validateFDboost(final)
fitted(final)
funMRD(final)
funMSE(final)
funRsquared(final)
# opt No.
folds <- cv(weights = rep(1, final$ydim[1]),type="subsampling",B=10)
opt <- applyFolds(final,folds=folds,grid = 1:100,mc.cores=2)
final <- final[mstop(opt)]
final <- final[10]
## plot the effect of temperature
plot(final, which = 1,pers = TRUE, main = "intercept with Poisson model")
plot(final, which = 2,pers = TRUE, main = "temperature", zlab = "")
plot(final, which = 2,pers = F, main = "temperature", zlab = "",col=hcl.colors(20,"YlGnBu"))
plot(final, which = 3,pers = TRUE, main = "humidity", zlab = "")
plot(final, which = 3,pers = F, main = "humidity", zlab = "",col=hcl.colors(20,"YlGnBu"))
plot(final, which = 4,pers = TRUE, main = "work")
# residuals
residual <- residuals(final)

# predict - work ----------------------------------------------------------

pred0 <- data__all
pred0$count <- scale(data__all$count,scale = F)
predict <- FDboost(work ~ 1 
                   + bsignal(temp, hour.s, knots = 10, df = 3)
                   + bsignal(humidity, hour.s, knots = 10,  df = 3)
                   + bsignal(count, hour.s, knots = 10, df = 3), 
                   timeformula = NULL,  control=boost_control(mstop = 1000), 
                   data=pred0,
                   family=Binomial())
pred <- predict(predict, type="response")
round_predes <- round(pred)
table(round_predes, as.numeric(pred0$work))


# plotting pred,res,obs ---------------------------------------------------

# ind <- sapply(1:31, function(s){ which(data_my$day == ord[s]) })
ind <- rep(1:20,1)
smoothRes <- predict(final)
residual <- residuals(final)
# if( is.null(dim(smoothRes)) ) smoothRes <- matrix(0, ncol = 24, nrow = 31)
# smoothRes <- (smoothRes)[ind, ]
# smoothRes <- ( predict(mod4, which=3) )[ind, ]
workOrd <- data__all$work[ind]
fit3 <- (predict(final))[ind, ]
response <- data__all$count[ind, ]
date <- data__all$date
library(maps)
par(mfrow=c(4,5))
for(i in 1:20) {
  plot(1:24, smoothRes[i, ], col = as.numeric(workOrd[i]), type="b",
       main = paste(date[i]),
       cex = 1.2, cex.axis = .8, ylab = "", xlab = "")
  abline(h = 0, col = 8)
}
# residuals
par(mfrow=c(4,5))
#layout(rbind(matrix(1:28, 4, 8), rep(29, 4), rep(29, 4)))
for(i in 1:20) {
  plot(1:24, residual[i, ], col = as.numeric(workOrd[i]),
       # ylim = range(residual_Jan, response-fit3),
       main = paste(data__all$date[i]),
       cex = 1.2, cex.axis = .8, ylab = "", xlab = "")
  abline(h = 0, col = 8)
}
# observed
par(mfrow=c(4,5))
#layout(rbind(matrix(1:28, 4, 8), rep(29, 4), rep(29, 4)))
for(i in 1:20) {
  plot(1:24, data__all$count[i, ], col = as.numeric(workOrd[i]), type="b",
       # ylim = range(residual_Jan, response-fit3),
       main = paste(data__all$date[i]),
       cex = 1.2, cex.axis = .8, ylab = "", xlab = "")
  abline(h = 0, col = 8)
}
# plotPredicted(final_Jan)


