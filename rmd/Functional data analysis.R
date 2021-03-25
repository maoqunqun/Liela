# Set up ------------------------------------------------------------------
library(readr)
library(fdapace)
library(tidyverse)
library(dplyr)
library(plyr)
library(fda)
library(lubridate)
library(ggplot2)
library(fdapace)
rm(list=ls())
hour <- read.csv("hour.csv")
data0 <- hour %>% filter(yr==0&mnth==1) %>% 
  mutate(day=day(dteday)) %>%  
  dplyr::select(hr,cnt,day)  

data <- data0 %>% pivot_wider(
  names_from = day,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.data.frame()
# Smoothing - 1D local linear kernel---------------------------------------------------------------
## fit curves ##
# original curves
hour_df = hour %>% filter(yr == 0& mnth==1) %>% 
  mutate(day=lubridate::day(dteday)) %>% 
  dplyr::select(instant, cnt, day) %>% 
  ddply(.(day), mutate, hours = 1:length(cnt), 
        count_smooth = Lwls1D(bw = 2, 
                              kernel_type = 'epan', 
                              xin = hours, 
                              yin = cnt, 
                              xout = hours))
local_linear <- hour_df %>% 
  dplyr::select(day,count_smooth,hours) %>%  
  pivot_wider(
    names_from = day,
    values_from = count_smooth,
    values_fill = list(count_smooth = 0)
  ) %>% as.data.frame()%>% dplyr::select(-hours)

hour_df_der = hour %>% filter(yr == 0& mnth==1) %>% 
  mutate(day=lubridate::day(dteday)) %>% 
  dplyr::select(instant, cnt, day) %>%  
  ddply(.(day), mutate, hours = 1:length(cnt), 
        count_smooth=Lwls1D(bw = 4, 
                            kernel_type = 'epan', 
                            xin = hours,
                            yin = cnt,
                            xout = hours,
                            nder=1))

local_linear_der <- hour_df_der %>% 
  dplyr::select(day,count_smooth,hours) %>%  
  pivot_wider(
    names_from = day,
    values_from = count_smooth,
    values_fill = list(count_smooth = 0)
  ) %>% as.data.frame()%>% dplyr::select(-hours)
# workingday=1
hour_df_der_w1 = hour %>% 
  filter(yr == 0& mnth==1 & workingday==1) %>% 
  mutate(day=lubridate::day(dteday)) %>% 
  dplyr::select(instant, cnt, day) %>%  
  ddply(.(day), mutate, 
        hours = 1:length(cnt), 
        count_smooth=Lwls1D(bw = 4,
                            kernel_type = 'epan', 
                            xin = hours,  yin = cnt,  
                            xout = hours,
                            nder=1))

local_linear_der_w1 <- hour_df_der_w1 %>% 
  dplyr::select(day,count_smooth,hours) %>%  
  pivot_wider(
    names_from = day,
    values_from = count_smooth,
    values_fill = list(count_smooth = 0)
  ) %>% as.data.frame() %>% dplyr::select(-hours)
# working day=0
hour_df_der_w0 = hour %>% 
  filter(yr == 0& mnth==1 & workingday==0) %>% 
  mutate(day=lubridate::day(dteday)) %>% 
  dplyr::select(instant, cnt, day) %>%  
  ddply(.(day), mutate, hours = 1:length(cnt), 
        count_smooth=Lwls1D(bw = 4, 
                            kernel_type = 'epan', 
                            xin = hours, 
                            yin = cnt, 
                            xout = hours,
                            nder=1))

local_linear_der_w0 <- hour_df_der_w0 %>% 
  dplyr::select(day,count_smooth,hours) %>%  
  pivot_wider(
    names_from = day,
    values_from = count_smooth,
    values_fill = list(count_smooth = 0)
  ) %>% as.data.frame()%>% dplyr::select(-hours)
## ploting ##
quartz()
par(mfrow=c(2,2))
# orginal
matplot(0:23,local_linear, type = "l", lty = 1, xlab='hour', 
        ylab= 'counts (y)',main = list("(a): Smoothed original curves", 
                                       cex = 0.8, font = 4))
# one derivative
matplot(0:23,local_linear_der, type = "l", lty = 1, xlab='hour',
        ylab='y\'',main = list("(b): Smoothed first derivative curves", 
                               cex = 0.8, font = 4))
# working
matplot(0:23,local_linear_der_w1, type = "l", lty = 1, xlab='hour',
        ylab = '',main = list(
          "(c): Smoothed first derivative curves for working days", 
          cex = 0.8, font = 4))
# non-working
matplot(0:23,local_linear_der_w0, type = "l", lty = 1, xlab='hour',
        ylab = '',main = list(
          "(d): Smoothed first derivative curves for non-working days", 
          cex = 0.8, font = 4))
# Smoothing - Bspline -----------------------------------------------------
tobs = data[,1]
nobs = length(tobs)
knots    = c(seq(0,23,1));
nknots   = length(knots);
norder   = 4;
nbasis   = length(knots) + norder - 2;
basis    = create.bspline.basis(c(min(tobs),max(tobs)),
                                nbasis,norder,knots);
basismat = eval.basis(tobs, basis);
# Use quadrature to get integral - Composite Simpson's Rule
delta <- 0.02
quadpts <- seq(0,1,delta)
nquadpts <- length(quadpts)
quadwts <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
quadwts <- c(1,rep(c(4,2),(nquadpts-1)/2))
quadwts[nquadpts] <- 1
quadwts <- quadwts*delta/3
# Second derivative of basis functions at quadrature points
Q2basismat = eval.basis(quadpts, basis,2);
# estimates for basis coefficients
Rmat = t(Q2basismat)%*%(Q2basismat*(quadwts%*%t(rep(1,nbasis))))
# estimates for basis coefficients
basismat2 = t(basismat)%*%basismat;
lambda = 0.8   # smoothing parameter
Bmat  = basismat2 + lambda*Rmat;
data1<-data %>% dplyr::select(-hr) %>% as.matrix()
chat = ginv(Bmat)%*%t(basismat)%*%data1;
yhat = basismat%*%chat;
yhat2 = basismat%*%ginv(t(basismat)%*%basismat)%*%t(basismat)%*%data1 
quartz()
matplot(0:23, yhat, type = "l", lty = 1, xlab='hour', ylab= 'counts (y)',
        main = list("B-Spline smoothed original curves", cex = 0.8, font = 4))
# FPCA --------------------------------------------------------------------
a<- hour %>% filter(yr == 0& mnth==1) %>% mutate(day=lubridate::day(dteday)) %>% dplyr::select(hr, cnt, day)
abc <- MakeFPCAInputs(a$day, a$hr, a$cnt)
# mu=2, cov=2.5
fpca <- FPCA(abc$Ly,abc$Lt,list(plot = TRUE, userBwMu=2, userBwCov=2.5))
# atuo - mu=1.15, ocv=2.3
fpca_auto <- FPCA(abc$Ly,abc$Lt,list(plot = TRUE, methodMuCovEst="smooth"))
fpca_auto$bwMu;fpca_auto$bwCov
#compare
quartz()
par(mfrow=c(1,2))
CreatePathPlot( fpca_auto, subset = 1:10, main = "(a): auto smooth", pch = 16)
CreatePathPlot( fpca, subset = 1:10, main = " (b): mu=2, cov=2.5", pch = 16)
# fpca
quartz()
plot(fpca)
CreateOutliersPlot(fpca, optns = list(K = 3, variant = 'KDE'))
CreateFuncBoxPlot(fpca, xlab = 'Hours', ylab = '# Renting', optns = list(variant='bagplot'))
SelectK(fpca,FVEthreshold = 0.99)
#fpca - derivative
der0 <- fitted(fpca,derOptns = list(p=1)) 
der <- t(der0)%>%  as.data.frame() 
matplot(0:23,der, type = "l", lty = 1, main="FPCA: First derivation")
# FCReg -------------------------------------------------------------------
# data frame set up
data0 <- hour %>% filter(yr==0&mnth==1) %>% mutate(day=day(dteday)) %>%  
  dplyr::select(hr,cnt,day)  
data1 <- data0 %>% pivot_wider(
  names_from = day,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.matrix()
data2 <-hour %>% filter(yr==0&mnth==1) %>% mutate(day=day(dteday)) %>%  
  dplyr::select(hr,atemp,day)  %>% pivot_wider(
    names_from = day,
    values_from = atemp,
    values_fill = list(atemp = 0)
  ) %>% as.matrix()
# 31*24
counting <- t(data1[,2:32] )
tempture <- t(data2[,2:32])
t <- seq(0,23,by=1)
X_1 <- Sparsify(tempture,t,sparsity =24)
Y_1 <- Sparsify(counting,t,sparsity =24) #c(10:20)
vars_ <-list(X_1,Y=Y_1)
fcreg <- FCReg(vars_,2,2.5,t)
# visulization
intercept <- fcreg$beta0 %>% as.data.frame() %>% 
  rename_at( vars(starts_with(".")),~ str_replace(., ".", "beta0"))%>% 
  mutate(out=c(0:23))
slope <- t(fcreg$beta)%>% as.data.frame() %>% 
  rename_at( vars(starts_with("V")),~ str_replace(., "V1", "beta"))%>%
  mutate(out=c(0:23))
visual <- intercept %>% inner_join(slope)
col=c(rep(1,9),rep(2,8),rep(1,7))
data_mean <- data0 %>% dplyr::group_by(hr) %>% dplyr::summarize(mean=mean(cnt))
attach(visual)
quartz()
par(mfrow=c(1,3))
plot(out,beta0,col=col, main="(a): Plot of intercept by hour", xlab = "hour");lines(out,beta0,col="grey80")
plot(out,beta,col=col, main="(b): Plot of slope (temp) by hour", xlab = "hour");lines(out,beta,col="grey80")
plot(data_mean$hr,col=col,data_mean$mean,main="(c): Plot of mean count per hour for all days (all)", xlab= "hour", ylab="mean");lines(data_mean$hr,data_mean$mean,col="grey80")
detach(visual)


# working==0 --------------------------------------------------------------

# data frame set up
data0 <- hour %>% filter(yr==0&mnth==1&workingday==0) %>% mutate(day=day(dteday)) %>%  dplyr::select(hr,cnt,day)  
data1 <- data0 %>% pivot_wider(
  names_from = day,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.matrix()
data2 <-hour %>% filter(yr==0&mnth==1) %>% mutate(day=day(dteday)) %>%  dplyr::select(hr,atemp,day)  %>% pivot_wider(
  names_from = day,
  values_from = atemp,
  values_fill = list(atemp = 0)
) %>% as.matrix()

# 31*24
counting <- t(data1[,2:11] )
tempture <- t(data2[,2:11])
t <- seq(0,23,by=1)
X_1 <- Sparsify(tempture,t,sparsity =24)
Y_1 <- Sparsify(counting,t,sparsity =24) #c(10:20)
vars_ <-list(X_1,Y=Y_1)
fcreg <- FCReg(vars_,2,2.5,t)
# visulization
intercept <- fcreg$beta0 %>% as.data.frame() %>% rename_at( vars(starts_with(".")),  ~ str_replace(., ".", "beta0")) %>% mutate(out=c(0:23))
slope <- t(fcreg$beta)%>% as.data.frame() %>% rename_at( vars(starts_with("V")),  ~ str_replace(., "V1", "beta"))  %>%mutate(out=c(0:23))
visual <- intercept %>% inner_join(slope)
col=c(rep(1,9),rep(3,8),rep(1,7))
data_mean <- data0 %>% dplyr::group_by(hr) %>% dplyr::summarize(mean=mean(cnt))
attach(visual)
quartz()
par(mfrow=c(1,3))
plot(out,beta0,col=col, main="(a): Plot of intercept by hour", xlab = "hour");lines(out,beta0,col="grey80")
plot(out,beta,col=col, main="(b): Plot of slope (temp) by hour", xlab = "hour");lines(out,beta,col="grey80")
plot(data_mean$hr,col=col,data_mean$mean,main="(c): Plot of mean count per hour for all days (non-working)", xlab= "hour", ylab="mean");lines(data_mean$hr,data_mean$mean,col="grey80")
detach(visual)

# working==1 --------------------------------------------------------------

# data frame set up
data0 <- hour %>% filter(yr==0&mnth==1&workingday==1) %>% mutate(day=day(dteday)) %>%  dplyr::select(hr,cnt,day)  
data1 <- data0 %>% pivot_wider(
  names_from = day,
  values_from = cnt,
  values_fill = list(cnt = 0)
) %>% as.matrix()
data2 <-hour %>% filter(yr==0&mnth==1) %>% mutate(day=day(dteday)) %>%  dplyr::select(hr,atemp,day)  %>% pivot_wider(
  names_from = day,
  values_from = atemp,
  values_fill = list(atemp = 0)
) %>% as.matrix()

# 31*24
counting <- t(data1[,2:20] )
tempture <- t(data2[,2:20])
t <- seq(0,23,by=1)
X_1 <- Sparsify(tempture,t,sparsity =24)
Y_1 <- Sparsify(counting,t,sparsity =24) #c(10:20)
vars_ <-list(X_1,Y=Y_1)
fcreg <- FCReg(vars_,2,2.5,t)
# visulization
intercept <- fcreg$beta0 %>% as.data.frame() %>% rename_at( vars(starts_with(".")),  ~ str_replace(., ".", "beta0")) %>% mutate(out=c(0:23))
slope <- t(fcreg$beta)%>% as.data.frame() %>% rename_at( vars(starts_with("V")),  ~ str_replace(., "V1", "beta"))  %>%mutate(out=c(0:23))
visual <- intercept %>% inner_join(slope)
col=c(rep(1,9),rep(4,8),rep(1,7))
data_mean <- data0 %>% dplyr::group_by(hr) %>% dplyr::summarize(mean=mean(cnt))
attach(visual)
quartz()
par(mfrow=c(1,3))
plot(out,beta0,col=col, main="(a): Plot of intercept by hour ", xlab = "hour");lines(out,beta0,col="grey80")
plot(out,beta,col=col, main="(b): Plot of slope (temp) by hour", xlab = "hour");lines(out,beta,col="grey80")
plot(data_mean$hr,col=col,data_mean$mean,main="(c): Plot of mean count per hour for all days (working)", xlab= "hour", ylab="mean");lines(data_mean$hr,data_mean$mean,col="grey80")
detach(visual)


