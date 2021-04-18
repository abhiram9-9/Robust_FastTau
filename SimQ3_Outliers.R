#clean

rm(list = ls())

#-----------------------

source("FastTau/FastTau.R")

n = 100

#-----------


library(MASS)
library(latex2exp)

makeData = function(p, eps, m){
  
  # '.g' for good, '.b' for bad
  
  n.g = (1-eps)*n
  n.b = n - n.g
  
  mu.g = rep(0, p-1)
  mu.b = c(100, rep(0,p-2))
  I = diag(rep(1, p-1))
  
  #good data
  X.g = mvrnorm(n = n.g, mu = mu.g, Sigma = I)
  X.g.p = as.matrix(rep(1, n.g))          #pth covariate : x_ip = 1
  X.g = cbind(X.g, X.g.p) 
  
  err = rnorm(n = n.g, mean = 0, sd = 1)
  B_0 = as.matrix(rep(0, p))              #coefficient vector
  y.g =  X.g %*% (B_0) + err
  
  #bad data : high leverage points
  X.b = mvrnorm(n = n.b, mu = mu.b, Sigma = 0.1^2 * I)
  X.b = cbind(X.b, as.matrix(rep(1, n.b)))  #pth covariate : x_ip = 1
  
  y.b = as.matrix(rnorm(n = n.b, mean = m, sd = 0.1^2))
  
  
  #full dataset
  
  X = rbind(X.g, X.b)         # n X (p-1) matrix
  y = rbind(y.g, y.b)
  
  return(list("X" = X, "y" = y))
}  

# p = {2, 5, 10}
# eps = {0.1, 0.2}
# m = {100, 120, ... 200}



########################################
p = 2
eps = 0.2
m = 120

data = makeData(p, eps, m)
dim(data$X)
dim(data$y)
head(data$X)
X = data$X
y = data$y


#bringing the column of 1's to the front
X = X[,c(p, 1:p-1)]

pairs(X[,-1], pch = 19)        #invalid for p=2

#X without the 1's column
X.No1 = as.matrix(X[,-1])

#exploring the dataset
library(rrcov)
library(robustbase)
library(corrplot)

m = colMeans(X.No1)
S = cov(X.No1)

col3 <- colorRampPalette(c("red", "white", "blue"))
corrplot(S, is.corr = FALSE, method = "color", col = col3(10), 
         tl.col= "black", type = 'lower', tl.cex = 1.2, cl.lim = c(min(S), max(S)),
         addgrid.col = "grey", tl.srt=0)

X.mcd = CovMcd(X.No1) 
plot(X.mcd,which = "tolEllipse",classic=TRUE)
covPlot(X.No1,which="dd",m.cov=X.mcd)

plot(X.mcd, which = "distance", classic = TRUE)# 2 plots
plot(X.mcd, which = "dd")

#for p = 2
plot(X.No1, y)
abline(h = mad(X.No1), col = "red", lwd = 2)


library(mrfDepth)
fit.s.d <- outlyingness(X.No1)
which(fit.s.d$outlyingnessX >= fit.s.d$cutoff)
plot(fit.s.d$outlyingnessX, pch = 20, cex = 1.3, ylab = "Outlyingness", xlab = "Index",
     cex.lab = 1.3, cex.axis = 1.3); grid()
abline(h = fit.s.d$cutoff, lwd = 3, col  ="red")

D2 <- sqrt(mahalanobis(X.No1, m, S))            #Claasical Mah. Dist.
plot(D2, pch = 19)
abline(h = sqrt(qchisq(p = 0.975, df = 1)), col = "red", lwd = 2)

plot(D2, fit.s.d$outlyingnessX, ylab = "Robust(SDO) distances", xlab = "Mah. distances", cex = 1.3)


#ESTIMATION

#-------------------------------------------------

lmOutPlot <- function(standardResid, dat, datMean, datCov){
  
  MDLS <- sqrt(mahalanobis(dat, datMean, datCov))
  plot(MDLS,standardResid, xlab="MD", pch = 19, ylab="standardised residual", main="standardised residuals vs MD",ylim=c(-10,10) )
  abline(h = 2.5, col = "red")
  abline(h = -2.5, col = "red")
  abline(v= sqrt(qchisq(0.975,ncol(dat))), col = "red")
  identify(MDLS, standardResid,plot=T)
  
}


###########
#  LS     #
###########


LSModel <- lm(y ~ X.No1, x = T)
LSModel$coefficients
standardlm <- rstandard(LSModel)
#meanAlc <- colMeans(X.No1)
#covAlc <- cov(X.No1)

lmOutPlot(standardlm,X.No1, m, S)         #classical mah. distance
lmOutPlot(standardlm,X.No1, X.mcd$center, X.mcd$cov)         #Robust(CovMCD) mah. distance

names(summary(LSModel))
names(LSModel)

summary(LSModel)$sigma

###########
#  LTS     #
###########
?ltsReg
LTSModel <- ltsReg(y ~ X.No1, mcd=T, x.ret = T) #, mcd=TRUE) #intercept automatically added

lmOutPlot(LTSModel$residuals/LTSModel$scale,X.No1, X.mcd$center, X.mcd$cov)

#plot(LTSModel)
plot(LTSModel,which="rdiag")
names(summary(LTSModel))
names(LTSModel)
LTSModel$scale
summary(LTSModel)$sigma

###########
#  M     #
###########

?lmrob..M..fit
m1 <- lmrob..M..fit(X, y, beta.initial = coef(LTSModel), scale = LTSModel$scale,
                    control = lmrob.control(copute.rd=T)) #dataset should have an intercept column
m1. = rlm(X, y, method = "M")
names((m1.))
m1.$scale

m1$coefficients

names(m1$control)
m1$control$tuning.chi
m1$control$method

lmOutPlot(m1$residuals/m1[["scale"]],X.No1, X.mcd$center, X.mcd$cov)

max(abs(m1$residuals/m1[["scale"]]))## increase y lim in outliermap
names(m1)
m1$scale

###########
#  S     #
###########

Ops <- lmrob.control(copute.rd=T)
Ops$max.it <- 5000
SModel <- lmrob(y ~ X.No1, control = Ops)[["init.S"]]

lmOutPlot(SModel$residuals/SModel[["scale"]],X.No1,X.mcd$center, X.mcd$cov)
names(SModel)
SModel$scale

###########
#  MM     #
###########
Ops <- lmrob.control(copute.rd=T)
Ops$max.it <- 10000
MMModel <- lmrob(y ~ X.No1, control = Ops)
plot(MMModel,which=1)
lmOutPlot(MMModel$residuals/MMModel$scale, X.No1,X.mcd$center, X.mcd$cov)


names(MMModel)
MMModel$scale
names(summary(MMModel))
summary(MMModel)$sigma

#####################################################
# Regression of median value by t est Tau estimator #
#####################################################
ptm <- proc.time()
tauest <- FastTau(x=X, y=y, N=500, kk=2, tt=5, rr=2, approximate=0, seed=456)
lmOutPlot(MMModel$residuals/MMModel$scale, X.No1,X.mcd$center, X.mcd$cov)

proc.time() - ptm
names(tauest)
tauest$beta

###########
#  coefs  #
###########

coefsMat <- cbind(LSModel$coefficients,LTSModel$coefficients, m1$coefficients, SModel$coefficients,MMModel$coefficients, tauest$beta)
colnames(coefsMat) <- c("LS", "LTS", "M", "S", "MM", "Tau")


corrplot(coefsMat, is.corr = FALSE, method = "color",  
         tl.col= "black", tl.cex = 1.2,
         addgrid.col = "grey")


#################
#  scale/sigma  #
#################


scale.list = as.matrix(c(summary(LSModel)$sigma, LTSModel$scale, m1$scale, SModel$scale, 
                         MMModel$scale, tauest$scale))
rownames(scale.list) = c("LS", "LTS", "M", "S", "MM", "Tau")
colnames(scale.list) = c("Scale")


